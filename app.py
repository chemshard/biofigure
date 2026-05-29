from flask import Flask, render_template, request, jsonify, Response
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote, unquote
from bs4 import BeautifulSoup
import requests
import threading
import time
import re
import random

app = Flask(__name__)

# ─────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────

IBO_TOPICS = [
    "cell biology", "genetics", "molecular biology",
    "photosynthesis", "cellular respiration", "mitosis meiosis",
    "DNA replication", "protein synthesis", "evolution", "ecology",
    "animal physiology", "plant physiology", "neuroscience",
    "immunology", "microbiology", "biochemistry enzymes",
    "membrane transport", "gene expression", "hormone signaling",
    "osmoregulation"
]

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0 Safari/537.36"
    )
}

# ── NCBI rate-limit guard ──────────────────────────────────────
# Without an API key: 3 req/s max.  With an API key: 10 req/s.
# We use a semaphore + inter-request delay to stay well under.
# Set NCBI_API_KEY = "" to use the anonymous tier (slower).
NCBI_API_KEY = ""          # paste your key here if you have one
_NCBI_RPS    = 2.5 if not NCBI_API_KEY else 9   # requests per second
_ncbi_sem    = threading.Semaphore(3 if not NCBI_API_KEY else 8)
_last_ncbi   = [0.0]
_ncbi_lock   = threading.Lock()

def ncbi_get(url, params=None, **kwargs):
    """
    Throttled GET to any NCBI Entrez endpoint.
    - Enforces a minimum inter-request interval.
    - Retries up to 3× on 429 with exponential back-off.
    - Automatically appends api_key when configured.
    """
    if NCBI_API_KEY:
        params = dict(params or {})
        params["api_key"] = NCBI_API_KEY

    kwargs.setdefault("headers", HEADERS)
    kwargs.setdefault("timeout", 20)

    for attempt in range(4):
        with _ncbi_sem:
            # Enforce minimum gap between requests
            with _ncbi_lock:
                gap = 1.0 / _NCBI_RPS
                elapsed = time.time() - _last_ncbi[0]
                if elapsed < gap:
                    time.sleep(gap - elapsed)
                _last_ncbi[0] = time.time()

            try:
                r = requests.get(url, params=params, **kwargs)
            except requests.RequestException as e:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
                continue

        if r.status_code == 429:
            wait = 2 ** (attempt + 1)
            print(f"[NCBI 429] backing off {wait}s (attempt {attempt+1})")
            time.sleep(wait)
            continue
        return r

    return r   # return whatever we have after all retries

_img_cache = {}

# ─────────────────────────────────────────────────────────────
# Step 1: Search PMC via Entrez
# ─────────────────────────────────────────────────────────────

def search_pmc(query, max_results=8):
    try:
        r = ncbi_get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                "db":      "pmc",
                "term":    f"{query}[Title/Abstract] AND open access[filter]",
                "retmax":  50,          # fetch a large pool to randomise from
                "retmode": "json",
                "sort":    "relevance", # keep relevance so pool stays on-topic
            },
        )
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        print(f"[esearch error] {e}")
        return []

    if not ids:
        return []

    # Shuffle so each search surfaces different articles from the relevant pool
    random.shuffle(ids)
    ids = ids[: max_results * 3]   # trim after shuffle

    try:
        r = ncbi_get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "pmc", "id": ",".join(ids), "retmode": "json"},
        )
        r.raise_for_status()
        summary = r.json().get("result", {})
    except Exception as e:
        print(f"[esummary error] {e}")
        summary = {}

    papers = []
    for uid in ids:
        doc = summary.get(uid, {})
        if not doc or uid == "uids":
            continue
        pmc_id  = doc.get("uid") or uid
        title   = doc.get("title", f"PMC{pmc_id}").rstrip(".")
        authors = ", ".join(a.get("name", "") for a in doc.get("authors", [])[:3])
        if len(doc.get("authors", [])) > 3:
            authors += " et al."
        journal = doc.get("fulljournalname") or doc.get("source", "")
        year    = (doc.get("pubdate") or "")[:4]
        papers.append({
            "pmc_id":  pmc_id,
            "title":   title,
            "authors": authors,
            "journal": journal,
            "year":    year,
            "url":     f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/",
        })
    return papers

# ─────────────────────────────────────────────────────────────
# Step 2a: Scrape blob URLs from the PMC article's inline JSON
#
# PMC embeds a JSON block in the static HTML (no JS needed):
#   <script id="article-page-data" type="application/json">
#     { ..., "figures": [ { "id": "F1", "src": "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/..." }, ... ] }
#   </script>
#
# We also sweep every string in that JSON for any CDN blob URL
# and key each one by its filename stem so we can match it to
# the JATS <graphic xlink:href> value from eFetch.
#
# Real URL structure (from observed URLs):
#   https://cdn.ncbi.nlm.nih.gov/pmc/blobs/{h1}/{pmc_numeric}/{h2}/{stem}.jpg
# ─────────────────────────────────────────────────────────────

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Matches full (https://) and protocol-relative (//cdn...) blob URLs.
# Hex segment lengths are relaxed — PMC CDN path structure varies by article vintage.
_BLOB_RE   = re.compile(
    r'(?:https?:)?//cdn\.ncbi\.nlm\.nih\.gov/pmc/blobs/'
    r'[0-9a-f]+/\d+/[0-9a-f]+/[^\s"\'<>&?#]+',
    re.I,
)


def _normalize_blob_url(url):
    """Ensure blob URL has an explicit https: scheme."""
    return url if url.startswith("http") else "https:" + url


def _strip_img_ext(s):
    """Remove a trailing image extension from a stem string (handles .jpg.jpg too)."""
    return re.sub(r'(\.[a-z]{2,4})+$', '', s, flags=re.I)


def _scan_html_for_blobs(html):
    """
    Scan a page of HTML for CDN blob URLs. Returns {stem: url}.

    Pass 1: <script type="application/json"> blocks — fastest, most reliable.
    Pass 2: all other <script> tags (plain JS / webpack bundles).
    Pass 3: <img> and <source> attributes (src, data-src, data-lazy-src, srcset).
    Pass 4: full-page regex sweep — last resort.
    """
    blob_map = {}

    def _add(url):
        blob_url = _normalize_blob_url(url)
        stem = _strip_img_ext(blob_url.rstrip("/").split("/")[-1])
        if stem and stem not in blob_map:
            blob_map[stem] = blob_url

    soup = BeautifulSoup(html, "lxml")

    # Pass 1: typed JSON script blocks
    for script in soup.find_all("script", type="application/json"):
        for m in _BLOB_RE.finditer(script.string or ""):
            _add(m.group(0))

    # Pass 2: untyped / text/javascript script blocks
    if not blob_map:
        for script in soup.find_all("script"):
            if script.get("type", "").lower() in ("", "text/javascript"):
                for m in _BLOB_RE.finditer(script.string or ""):
                    _add(m.group(0))

    # Pass 3: img / source element attributes
    if not blob_map:
        for tag in soup.find_all(["img", "source"]):
            for attr in ("src", "data-src", "data-lazy-src", "srcset"):
                for m in _BLOB_RE.finditer(tag.get(attr, "")):
                    _add(m.group(0))

    # Pass 4: full-page regex sweep
    if not blob_map:
        for m in _BLOB_RE.finditer(html):
            _add(m.group(0))

    return blob_map


def scrape_blob_urls(pmc_id):
    """
    Fetch the PMC article page (static HTML) and return a dict of
    { stem -> blob_url }.

    Tries two pages in order, stopping as soon as any blob URLs are found:
      1. Main article page  -- new-renderer articles embed all blob URLs here.
      2. /figures/ gallery  -- server-rendered gallery page; reliably contains
                               static blob URLs for articles whose main page
                               loads images entirely via JavaScript.

    stem = filename without any extension, e.g. "tpae135f1"
    """
    pages_to_try = [
        f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/",
        f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/figures/",
    ]

    for page_url in pages_to_try:
        try:
            r = requests.get(page_url, headers=HEADERS, timeout=15, allow_redirects=True)
            if r.status_code != 200:
                print(f"[PMC{pmc_id}] HTTP {r.status_code} for {page_url}")
                continue
        except Exception as e:
            print(f"[PMC{pmc_id}] fetch error ({page_url}): {e}")
            continue

        blob_map = _scan_html_for_blobs(r.text)
        if blob_map:
            print(f"[PMC{pmc_id}] blob_map: {len(blob_map)} URLs found via {page_url}")
            return blob_map

    print(f"[PMC{pmc_id}] blob_map: 0 URLs found")
    return {}


# ─────────────────────────────────────────────────────────────
# Step 2b: Get figure metadata from eFetch JATS XML
#
# Returns list of dicts: {fig_id, stem, label, caption}
#   fig_id  — the <fig id="..."> attribute, used to build the figure subpage URL
#   stem    — graphic href stem, used to match against blob_map keys
# ─────────────────────────────────────────────────────────────

def fetch_figures_from_xml(pmc_id):
    """Return list of {fig_id, stem, direct_url, label, caption} from eFetch XML.

    direct_url is populated when the graphic href is already a full CDN blob URL.
    Many articles embed it directly in the XML — no scraping or subpage needed.
    stem is always set for blob_map matching (PATH A).
    """
    try:
        r = ncbi_get(
            EFETCH_URL,
            params={"db": "pmc", "id": pmc_id, "rettype": "xml"},
        )
        if r.status_code != 200:
            print(f"[PMC{pmc_id}] eFetch HTTP {r.status_code}")
            return []
    except Exception as e:
        print(f"[PMC{pmc_id}] eFetch error: {e}")
        return []

    try:
        soup = BeautifulSoup(r.content, "lxml-xml")
    except Exception:
        return []

    results = []
    for fig in soup.find_all("fig"):
        if fig.find_parent("fig"):
            continue

        fig_id = fig.get("id", "")   # e.g. "F1", "fig1", "RSOB180246F1"

        label_tag = fig.find("label")
        label = label_tag.get_text(strip=True) if label_tag else ""

        cap_tag = fig.find("caption")
        caption = re.sub(r"\s+", " ", cap_tag.get_text(" ", strip=True)) if cap_tag else ""

        graphic = fig.find("graphic")
        if not graphic:
            continue

        href = (
            graphic.get("xlink:href")
            or graphic.get("{http://www.w3.org/1999/xlink}href")
            or graphic.get("href")
            or ""
        ).strip()

        if not href:
            continue

        # Preserve full URL when the href is already a CDN blob URL.
        # Many articles (especially those whose HTML has no inline JSON) embed
        # the complete blob URL directly in the eFetch XML graphic element.
        direct_url = None
        if href.startswith("http"):
            direct_url = href
            filename   = href.rstrip("/").split("/")[-1]
            stem       = _strip_img_ext(filename)
        else:
            stem = _strip_img_ext(href.split("/")[-1])

        if stem:
            results.append({
                "fig_id":     fig_id,
                "stem":       stem,
                "direct_url": direct_url,
                "label":      label,
                "caption":    caption,
            })

    return results


# ─────────────────────────────────────────────────────────────
# Step 2c: Fetch blob URL from a figure's dedicated subpage
#
# For JS-rendered PMC articles, the main article HTML contains no blob URLs.
# Each figure has a server-rendered subpage that DOES contain them.
#
# PMC uses two different identifiers as the URL path component depending on
# the article renderer and vintage:
#   • fig_id  — the <fig id="..."> attribute, e.g. "F1", "fig1", "koab234-F1"
#   • stem    — the graphic href stem, e.g. "nihms-1821239-f0001"
#
# We try both (plus case variants and singular/plural path) so we cover all
# renderer combinations.
# ─────────────────────────────────────────────────────────────

def _blob_from_figure_subpage(pmc_id, fig_id, stem=None):
    """
    Fetch the per-figure subpage and return the first blob URL found, or None.

    Candidate URLs are built from both fig_id and stem (when provided) so that
    articles where the subpage key is the graphic href stem rather than the
    XML fig id are also resolved.
    """
    ids_to_try = []
    if fig_id:
        ids_to_try.append(fig_id)
    if stem and stem != fig_id:
        ids_to_try.append(stem)

    if not ids_to_try:
        return None

    seen = set()
    candidates = []
    for fid in ids_to_try:
        for fid_variant in (fid, fid.lower()):
            for path in ("figure", "figures"):
                url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/{path}/{fid_variant}/"
                if url not in seen:
                    seen.add(url)
                    candidates.append(url)

    for url in candidates:
        try:
            r = requests.get(url, headers=HEADERS, timeout=10, allow_redirects=True)
            if r.status_code != 200:
                continue
            blob_map = _scan_html_for_blobs(r.text)
            if blob_map:
                return next(iter(blob_map.values()))
        except Exception:
            continue
    return None


# ─────────────────────────────────────────────────────────────
# Step 2 (combined): blob URLs + figure metadata → figure list
#
# Three paths, tried in priority order:
#
# PATH A — HTML inline JSON has blob URLs (new renderer, post ~2022):
#   scrape_blob_urls() finds them in <script type="application/json">.
#   Match blob_map stems against eFetch XML figure stems.
#
# PATH B — eFetch XML graphic href is already a full CDN URL:
#   Many articles embed the complete blob URL directly in the XML.
#   We preserved it as direct_url — use it immediately, no scraping needed.
#
# PATH C — Last resort: fetch each figure's dedicated subpage:
#   /articles/PMC{id}/figure/{fig_id}/ is server-rendered and has the blob URL.
#   Only used when both blob_map is empty and direct_url is None.
# ─────────────────────────────────────────────────────────────

def get_figures_for_paper(paper, max_figs=5):
    pmc_id = paper["pmc_id"]
    jy = f"{paper['journal']} ({paper['year']})" if paper["journal"] else paper["year"]

    # Run HTML scrape and eFetch XML in parallel
    with ThreadPoolExecutor(max_workers=2) as ex:
        f_blobs  = ex.submit(scrape_blob_urls, pmc_id)
        f_xml    = ex.submit(fetch_figures_from_xml, pmc_id)
        blob_map = f_blobs.result()
        xml_figs = f_xml.result()

    if not xml_figs:
        print(f"[PMC{pmc_id}] skipping — no figures in XML")
        return []

    # Helper: build img-proxy URL with article-specific Referer hint.
    # The PMC CDN checks the Referer header; using the exact article page
    # instead of the bare domain root fixes 403s on blob URLs that were
    # discovered via figure subpages (Path C) rather than the main article page.
    def _proxy_url(img_url):
        article_ref = quote(f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmc_id}/", safe='')
        return f"/img-proxy?url={quote(img_url, safe='')}&ref={article_ref}"

    figures = []

    if blob_map:
        # ── PATH A: HTML inline JSON → match stems ────────────────
        path = "A"
        for fig in xml_figs:
            stem    = fig["stem"]
            img_url = blob_map.get(stem)

            if not img_url:
                for k, v in blob_map.items():
                    if k.startswith(stem) or stem.startswith(k):
                        img_url = v
                        break

            if not img_url:
                continue

            figures.append({
                "label":         fig["label"],
                "caption":       fig["caption"],
                "img_url":       _proxy_url(img_url),
                "pmc_id":        pmc_id,
                "paper_title":   paper["title"],
                "paper_url":     paper["url"],
                "paper_journal": jy,
            })
            if len(figures) >= max_figs:
                break

    elif any(f["direct_url"] for f in xml_figs):
        # ── PATH B: full CDN URL already in eFetch XML ────────────
        path = "B"
        direct_count = sum(1 for f in xml_figs if f["direct_url"])
        print(f"[PMC{pmc_id}] PATH B: {direct_count} direct URLs in XML")
        for fig in xml_figs:
            img_url = fig["direct_url"]
            if not img_url:
                continue
            figures.append({
                "label":         fig["label"],
                "caption":       fig["caption"],
                "img_url":       _proxy_url(img_url),
                "pmc_id":        pmc_id,
                "paper_title":   paper["title"],
                "paper_url":     paper["url"],
                "paper_journal": jy,
            })
            if len(figures) >= max_figs:
                break

    else:
        # ── PATH C: fetch each figure's dedicated subpage ─────────
        # /articles/PMC{id}/figure/{fig_id}/ is server-rendered and
        # contains the blob URL in static HTML.
        path = "C"
        target_figs = xml_figs[:max_figs]

        def _fetch_one(fig):
            return fig, _blob_from_figure_subpage(pmc_id, fig["fig_id"], fig["stem"])

        with ThreadPoolExecutor(max_workers=5) as ex:
            results = list(ex.map(_fetch_one, target_figs))

        for fig, blob_url in results:
            if not blob_url:
                print(f"[PMC{pmc_id}] no blob URL for fig_id={fig['fig_id']!r}")
                continue
            figures.append({
                "label":         fig["label"],
                "caption":       fig["caption"],
                "img_url":       _proxy_url(blob_url),
                "pmc_id":        pmc_id,
                "paper_title":   paper["title"],
                "paper_url":     paper["url"],
                "paper_journal": jy,
            })

    print(f"[PMC{pmc_id}] PATH {path} → {len(figures)} figures")
    return figures

# ─────────────────────────────────────────────────────────────
# Step 3: Orchestrate
# ─────────────────────────────────────────────────────────────

def fetch_figures_for_query(query, max_papers=5, max_figs_per_paper=5):
    papers = search_pmc(query, max_results=max_papers * 3)
    if not papers:
        return [], []

    all_papers  = []
    all_figures = []

    # Process papers in small batches. Each batch runs concurrently, but we
    # don't start the next batch until we know whether we still need more papers.
    # This avoids firing off 20+ fetches when 5 succeed in the first batch.
    batch_size = 3
    i = 0
    while i < len(papers) and len(all_papers) < max_papers:
        batch = papers[i : i + batch_size]
        i += batch_size

        with ThreadPoolExecutor(max_workers=batch_size) as ex:
            futures = {ex.submit(get_figures_for_paper, p, max_figs_per_paper): p for p in batch}
            for future in as_completed(futures):
                try:
                    figs = future.result()
                except Exception as e:
                    print(f"[future error] {e}")
                    continue
                if not figs:
                    continue
                all_papers.append(futures[future])
                all_figures.extend(figs)
                if len(all_papers) >= max_papers:
                    for f in futures:
                        f.cancel()
                    break

    print(f"[done] {len(all_papers)} papers, {len(all_figures)} figures total")
    return all_papers, all_figures

# ─────────────────────────────────────────────────────────────
# Image Proxy
# ─────────────────────────────────────────────────────────────

@app.route("/img-proxy")
def img_proxy():
    raw = request.args.get("url", "").strip()
    if not raw:
        return "", 400
    target = unquote(raw)
    if not target.startswith("http"):
        return "", 403

    # Use the article-specific Referer when provided; fall back to domain root.
    # The PMC CDN validates Referer and returns 403 if it doesn't match the
    # article that owns the image — especially for blob URLs found via subpages.
    raw_ref = request.args.get("ref", "").strip()
    referer = unquote(raw_ref) if raw_ref else "https://pmc.ncbi.nlm.nih.gov/"

    if target in _img_cache:
        data, ct = _img_cache[target]
        return Response(data, content_type=ct,
                        headers={"Cache-Control": "public, max-age=86400"})
    try:
        r = requests.get(
            target,
            headers={**HEADERS, "Referer": referer},
            timeout=20, stream=True, allow_redirects=True,
        )
        if r.status_code != 200:
            print(f"[img-proxy] {r.status_code} for {target}")
            return "", 404
        ct = r.headers.get("Content-Type", "")
        # Accept image/* and also application/octet-stream for some CDN blobs
        if not (ct.startswith("image/") or ct == "application/octet-stream"):
            # Try to guess from URL
            if not re.search(r'\.(jpg|jpeg|png|gif|webp|svg)$', target, re.I):
                return "", 404
        data = r.content
        if not ct.startswith("image/"):
            ct = "image/jpeg"   # safe default for CDN blobs
        _img_cache[target] = (data, ct)
        return Response(data, content_type=ct,
                        headers={"Cache-Control": "public, max-age=86400"})
    except Exception as e:
        print(f"[proxy error] {type(e).__name__}: {e}")
        return "", 502

# ─────────────────────────────────────────────────────────────
# Routes
# ─────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html", topics=IBO_TOPICS)

@app.route("/search")
def search():
    query = request.args.get("q", "").strip()
    if not query:
        return jsonify({"error": "No query provided"}), 400
    papers, figures = fetch_figures_for_query(query)
    return jsonify({
        "query":         query,
        "papers":        papers,
        "figures":       figures,
        "total_figures": len(figures),
        "total_papers":  len(papers),
    })

# ─────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app.run(debug=True, threaded=True, port=5000, use_reloader=False)