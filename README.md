# BioFigure

A Flask web app that pulls biology research paper figures and captions from
**Europe PMC** open-access articles. No API key required, completely free.

## Features
- Searches Europe PMC (covers all of PubMed/PMC) for open-access papers
- Uses the Europe PMC Figures API for structured figure + caption data
- Falls back to scraping PMC article pages if needed
- 20 pre-loaded IBO-relevant topic chips
- Click any figure to zoom in with full caption
- Links back to source papers

## Setup

```bash
# 1. Install dependencies (Python 3.8+)
pip install -r requirements.txt

# 2. Run
python app.py

# 3. Open in browser
http://localhost:5000
```

## How it works

1. **Search** — queries `https://www.ebi.ac.uk/europepmc/webservices/rest/search`
   with filters for open-access, full-text available papers
2. **Figures** — for each PMC article, fetches structured figure data from
   `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{id}/figures/json`
3. **Fallback** — if the figures API returns nothing, scrapes the PMC article
   HTML for `<figure>` elements and image URLs
4. **Display** — shows figures in a searchable grid with captions and paper links

## No API key needed
Uses only free, public APIs:
- Europe PMC REST API (ebi.ac.uk)
- PMC Open Access full-text (ncbi.nlm.nih.gov)
