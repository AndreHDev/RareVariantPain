import os

from .Caching import txt_cache
import requests
import xml.etree.ElementTree as ET
from time import sleep

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Rate limits (seconds between requests)
RATE_LIMIT_WITHOUT_KEY = 0.5
RATE_LIMIT_WITH_KEY = 0.2

API_KEY_FILE = os.path.join(os.path.dirname(__file__), "..", "API_key.txt")

def _load_api_key_from_file(path=API_KEY_FILE):
    try:
        with open(path, encoding="utf-8") as f:
            key = f.read().strip()
        return key or None
    except FileNotFoundError:
        return None

API_KEY = _load_api_key_from_file()

def _get_rate_limit(api_key=None):
    """Return the appropriate rate-limit delay (seconds) based on whether an API key is used."""
    return RATE_LIMIT_WITH_KEY if api_key else RATE_LIMIT_WITHOUT_KEY


def _request_with_retries(method, url, api_key=None, max_attempts=20, **kwargs):
    """Make an HTTP request with retry/backoff for transient errors.

    - Retries on network errors and 5xx server errors.
    - Uses increasing backoff (1s..5s) between retries.
    - Always waits the configured rate-limit delay before each attempt.
    """
    api_key = api_key or API_KEY

    for attempt in range(1, max_attempts + 1):
        sleep(_get_rate_limit(api_key))  # rate limit

        try:
            response = requests.request(method, url, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.HTTPError as e:
            status = getattr(e.response, "status_code", None)
            if status and status >= 500 and attempt < max_attempts:
                sleep(min(attempt, 5))
                continue
            raise
        except requests.exceptions.RequestException:
            if attempt < max_attempts:
                sleep(min(attempt, 5))
                continue
            raise


# Search pubmed with a general querry
@txt_cache()
def search_pubmed(query, retmax=500, api_key=None):
    api_key = api_key or API_KEY

    url = BASE_URL
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax
    }
    if api_key:
        params["api_key"] = api_key

    r = _request_with_retries("post", url, api_key=api_key, data=params)
    data = r.json()

    pmids_from_result = data["esearchresult"]["idlist"]
    if len(pmids_from_result) >= retmax:
        raise Exception(f"Result exceed maximum possible returnable results of {retmax}")

    return pmids_from_result

# Papers given paper cites
def get_references_from_pmid(pmid, api_key=None):
    api_key = api_key or API_KEY

    params = {
        "dbfrom": "pubmed",
        "db": "pubmed",
        "linkname": "pubmed_pubmed_refs",  
        "id": pmid,
        "retmode": "json"
    }
    if api_key:
        params["api_key"] = api_key

    r = _request_with_retries("get", ELINK_URL, api_key=api_key, params=params)
    data = r.json()
    linksets = data.get("linksets", [])
    if not linksets:
        return set()

    linksetdbs = linksets[0].get("linksetdbs", [])
    if not linksetdbs:
        return set()

    return set(linksetdbs[0].get("links", []))

# Papers that cite given paper
@txt_cache()
def get_citing_pmids(pmid, api_key=None):
    api_key = api_key or API_KEY

    params = {
        "dbfrom": "pubmed",
        "db": "pubmed", 
        "linkname": "pubmed_pubmed_citedin",
        "id": pmid,
        "retmode": "xml"
    }
    
    if api_key:
        params["api_key"] = api_key  

    response = _request_with_retries("get", ELINK_URL, api_key=api_key, params=params)

    root = ET.fromstring(response.text)

    citing_pmids = [
        link_id.text
        for link_id in root.findall(".//LinkSetDb[DbTo='pubmed']/Link/Id")
    ]

    return citing_pmids

# Get titles and abstracts for a list of PubMed IDs
@txt_cache()
def get_titles_and_abstracts(pmids, api_key=None):
    """
    Fetches titles and abstracts for a list of PubMed IDs.
    
    Args:
        pmids: List of PubMed IDs
        api_key: Optional NCBI API key for faster requests
        
    Returns:
        List of dictionaries with keys: pmid, title, abstract
    """
    api_key = api_key or API_KEY
    if not pmids:
        return []

    # Join pmids as comma-separated string
    id_string = ",".join(str(pmid) for pmid in pmids)
    
    params = {
        "db": "pubmed",
        "id": id_string,
        "retmode": "xml",
        "rettype": "medline"
    }
    
    if api_key:
        params["api_key"] = api_key
    
    response = _request_with_retries("get", EFETCH_URL, api_key=api_key, params=params)
    root = ET.fromstring(response.text)
    
    results = []
    for pubmed_article in root.findall(".//PubmedArticle"):
        # Extract PMID
        pmid_elem = pubmed_article.find(".//PMID")
        pmid = pmid_elem.text if pmid_elem is not None else None
        
        # Extract title
        title_elem = pubmed_article.find(".//ArticleTitle")
        title = title_elem.text if title_elem is not None else None
        
        # Extract abstract
        abstract_elem = pubmed_article.find(".//Abstract")
        if abstract_elem is not None:
            # Concatenate all abstract text elements
            abstract_parts = [elem.text for elem in abstract_elem.findall(".//AbstractText") if elem.text]
            abstract = " ".join(abstract_parts) if abstract_parts else None
        else:
            abstract = None
        
        results.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract
        })
    
    return results