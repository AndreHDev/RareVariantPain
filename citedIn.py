from APIs.PubMedAPI import *
import datetime
import os

# all scores from paper
score_pmids = [
    "38183205", "30371827", "24487276", "33618777", "25338716", "26727659",
    "26727659", "27666373", "25552646", "25552646", "36209109", "20354512",
    "12202775", "23056405", "21727090", "21152010", "15965027", "19858363",
    "25599402", "28288115", "24487584", "22955989", "22955989", "28968714",
    "25583119", "23033316", "28985712", "21457909", "30256891", "33893808",
    "24681721", "20676075", "28093075", "34717010", "12952881", "27193693",
    "33219223", "19734154", "15965030", "26301843", "30013180", "33462483",
    "16024819", "19478016", "17526529", "22241780", "22505138", "19602639",
    "16895930", "18387208", "16204125", "23819870", "26442818", "26213851",
    "29588361", "27197224", "11337480", "22689647", "12824425", "27776117",
    "4843792", "37733863", "1438297", "30661751", "30220433", "33479230",
    "30559491", "28498993", "27153718", "28851873", "27995669", "32352516",
    "34551312", "34861178", "37563329", "38934805", "35639618", "37636282",
    "40397311", "41286104", "34707284", "28092658", "34270679", "37550700",
    "37084271", "35453737", "32831124", "36376793", "29617928", "33300042",
    "37475887", "32340307", "35486646", "31096927", "30518757", "27569544",
    "33543123", "34349788", "27009626", "37198692", "35036922", "32735577",
    "33789710", "30823901", "30475984", "15285897", "14695534", "29689380",
    "35965322", "33761318", "37083939", "32101277", "30804562", "35079159",
    "36273432", "29750258", "34289339", "35801945", "36707993", "34648033",
    "36063453", "32938008", "26075791", "33686085", "32444882", "28592878",
    "36611253", "36651276", "31157880", "23736532", "28794409", "30704475",
    "34850938", "33866367", "33603233", "34608324", "39779956", "35551308",
    "35449021", "41606153", "35216679", "35817977", "26502339"
]
RVAT_pmids = [
    '22441326', '17101154', '18691683', '19810025', '19214210', 
    '20413981', '21072163', '21885029', '20471002', '21304886', 
    '21304886', '21737059', '19170135', '21408211', '22699862', 
    '23032573', '23483651', '23159251', '29932245', '29754769', 
    '21368279', '21368279', '18691683', '20413981', '29754769', 
    '30849328', '24831820', '32839606', '35653402', '21070896', 
    '22009789'
]

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