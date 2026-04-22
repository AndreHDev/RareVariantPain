from APIs.PubMedAPI import *
import datetime
import os

# all scores from paper
score_pmids = []
RVAT_pmids = []

# Pain–related query string
pain_query = """
( 
(chronic OR persisting OR persistent OR lasting OR neuropathic OR nociceptive OR nociplastic OR 
mixed OR neurogenic OR back OR neck OR migraine OR arthritis OR osteoart* OR joint OR rheumatic OR 
inflammatory OR musculoskeletal OR muscle OR visceral OR widespread OR somatoform OR cancer OR 
postoperative OR postsurgic* OR perioperative)
AND 
(pain OR painful OR orchialgia OR analgesi* OR fibromyalgia) 
)
"""


def get_citing_papers_matching_pain(pmids, query):
    # Create a timestamped results file inside the Results directory
    RESULTS_DIR = "Results"
    os.makedirs(RESULTS_DIR, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    results_fname = f"results_{ts}.txt"
    results_path = os.path.join(RESULTS_DIR, results_fname)
    print(f"Writing results to {results_path}")

    with open(results_path, "w", encoding="utf-8") as f:
        # Write initial info to the file
        print("Getting papers that cite given PMIDS and maching query.", file=f)
        print(f"Query: {query}", file=f)
        print(f"PMIDs to check: {pmids}", file=f)

        all_citing_matching_pain = {}
        citing_counts = []
        matching_counts = []

        for pmid in pmids:
            # For every PMID, get the papers that cite it
            print(f"Processing PMID {pmid}...")
            citing_pmids = get_citing_pmids(pmid)

            if not citing_pmids:
                print(f"PMID {pmid}: No citing papers found", file=f)
                citing_counts.append(0)
                matching_counts.append(0)
                continue

            # Build query to filter citing papers with pain query
            combined_query = query + " AND " + "(" + " ".join(citing_pmids) + ")"
            # Remove newlines for PubMed search
            combined_query = combined_query.replace("\n", "")  
            # Search PubMed with the combined query to find citing papers that match the pain query
            match_pmids = search_pubmed(combined_query)
            
            # Write results
            print(f"PMID {pmid}:", file=f)
            print(f"  Citing papers found: {len(citing_pmids)}", file=f)
            print(f"  Amount of citing papers matching pain query: {len(match_pmids)}", file=f)
            all_citing_matching_pain[pmid] = match_pmids
            citing_counts.append(len(citing_pmids))
            matching_counts.append(len(match_pmids))
            print("------", file=f)

        # Print summary lines
        print(f"Citing papers found: {' '.join(map(str, citing_counts))}", file=f)
        print(f"Matching pain query: {' '.join(map(str, matching_counts))}", file=f)

        # Collect all found pain PMIDs into a single set to remove duplicates. Makes manual review easier.
        all_pmids = set()
        for match_pmids in all_citing_matching_pain.values():
            all_pmids.update(match_pmids)  # add all PMIDs from the set
        # Join into a simple string with spaces
        all_pmids_str = " ".join(all_pmids)
        # Write to file
        print("All unique PMIDs from citing papers matching pain query:", file=f)
        f.write(all_pmids_str + "\n")

get_citing_papers_matching_pain(score_pmids, pain_query)