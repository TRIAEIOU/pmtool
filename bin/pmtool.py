import sys, argparse, json, os, re, math, requests, time, random
from bs4 import BeautifulSoup


def pmquery(q: str, rmax: int = -1):
    """Query PubMed"""
    if rmax < 0 or rmax > 100: qsize = 200
    elif rmax <= 10: qsize = 10
    elif rmax <= 20: qsize = 20
    elif rmax <= 50: qsize = 50
    elif rmax <= 100: qsize = 100
    qs = fr'https://pubmed.ncbi.nlm.nih.gov/?term={requests.utils.quote(q)}&size={qsize}&format=pubmed' # query string
    qn = 0 # query counter
    matchc = -1 # number of matches
    qc = -1 # query count to get recordc
    matches = ""
    rq = ""
    while True:
        r = requests.get(f"{qs}&page={qn + 1}")
        if r.status_code == 200:
            s = BeautifulSoup(r.text, 'html.parser')
            # First lap, calculate limits
            if qn == 0:
                matchc = int(s.select_one("meta[name=log_resultcount]")['content'])
                qc = math.ceil(min(rmax, matchc) / qsize)
                rq = s.select_one("meta[name='log_processedquery']")['content']

            recs = s.select_one('pre.search-results-chunk').text.replace("\r\n", "\n")
            # Last lap
            if qn == qc - 1:
                ok = False
                # Truncate if needed
                if rmax < matchc:
                    reca = recs.split('\n\n')
                    if rmod := rmax % qsize:
                        matches += "\n\n" + "\n\n".join(reca[:rmod])
                        ok = True
                if not ok:
                    matches += f"\n\n{recs}"
                break
            else:
                matches += f"\n\n{recs}"
                qn += 1
        else:
            break
    return {
        "user query": q,
        "actual query": rq,
        "result": matches.strip()
    }

def pmparse(input: str):
    """Parse string in PubMed format into dict"""
    articles = []
    for iarticle in input.split('\n\n'):
        iarticle = re.sub('[ ]\n[ ]+', ' ', iarticle)
        oarticle = {}
        authors = []
        for key, val in re.findall(r'^([A-Z]+)\s*- (.*?)$', iarticle, flags=re.MULTILINE):

            def add(dest, k, v):
                """Add key/val pair, if key exists convert val to list"""
                if not dest.get(k):
                    dest[k] = v
                elif type(dest[k]) != list:
                    dest[k] = [dest[k], val]
                else:
                    dest[k].append(val)

            if key == 'FAU':
                authors.append({'FAU': val})
            elif key == 'AU':
                authors[-1]['AU'] = val
            elif key == 'AUID':
                authors[-1]['AUID'] = val
            elif key == 'AD':
                if not authors[-1].get('AD'):
                    authors[-1]['AD'] = []
                authors[-1]['AD'].append(val)
            else:
                add(oarticle, key, val)

        oarticle['AUS'] = authors # Own key for array of authors
        oarticle['URL'] = f"https://pubmed.ncbi.nlm.nih.gov/{oarticle['PMID']}/" # Own key for URL
        articles.append(oarticle)
    
    return articles

def pmformat(article: dict, fmt: str = "md"):
    """Format article dict to string"""

    def stringify(val, sep = ", "):
        if type(val) == list:
            return sep.join(val)
        return val

    if a := article.get('AB'):
        abstract = "\n" + re.sub(
            r'(?:^| )([A-Z0-9 ]+:) ',
            lambda m: f"\n\n**{m.group(1).capitalize()}** ",
            a
        ).strip() + "\n"
    else:
        abstract = ""
    txts = []
    if article.get('PMC'):
        txts.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{article['PMC']}/")
    if lids := article.get('LID'):
        if type(lids) != list:
            lids = [lids]
        for lid in lids:
            parts = lid.rsplit(' ', 1)
            if len(parts) == 1:
                continue
            elif parts[1] == '[doi]':
                txts.append(f"https://doi.org/{parts[0]}")
            elif parts[1] == '[pii]':
                txts.append(f"https://linkinghub.elsevier.com/retrieve/pii/{parts[0]}")

    txt =  "  \n".join([f"[Full text]({_url})" for _url in txts])
    return f"""## {article['TI']}

Publication type: {stringify(article['PT'])}
{abstract}
### Information

Date of publication: {stringify(article['DP'])}  
Source: {stringify(article['SO'], ' | ')}  
Authors: {', '.join([au['AU'] for au in article['AUS']])}  
PMID: {article['PMID']}  
[PubMed entry]({article['URL']})  
{txt}"""

def main(argv):
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", "--input-file", type=str, help="read input from file (rather than stdin)")
    argparser.add_argument("-o", "--output-file", type=str, help="write output to file (rather than stdout)")
    argparser.add_argument("-f", "--format", type=str, help="output format (md/json)")
    argparser.add_argument("-n", "--number", type=int, help="max number of entries (mainly for use with queries)", default=-1)
    argparser.add_argument("-q", "--query", nargs=argparse.REMAINDER, help="interpret input as queries to run against PubMed")
    args = argparser.parse_args()

    input = ""
    if args.input_file:
        with open(args.input_file, encoding="utf-8") as f:
            for line in f:
                input += line
    else:
        for line in sys.stdin:
            input += line

    format = "md"
    if args.format:
        format = args.format
    elif args.output_file:
        if ext := os.path.splitext(args.output_file)[1]:
            ext = ext[1:]
            if ext in ['md', 'json']:
                format = ext

    if args.query != None:
        results = []
        for q in input.split('\n'):
            q = q.strip()
            if not len(q) or q[0] == '#':
                continue
            r = pmquery(q, args.number)
            r['result'] = pmparse(r["result"])
            results.append(r)
        
        if format == 'json':
            out = json.dumps(results)
        else:
            out = ""
            for r in results:
                out += f"# {r['user query']}\n\nActual query: {r['actual query']}\n\n"
                for rr in r['result']:
                    out += f"{pmformat(rr, 'md')}\n\n"
            out = out.strip()
    else:
        if format == 'json':
            out = json.dumps(pmparse(input))
        else:
            out = ""
            for r in pmparse(input):
                out += f"{pmformat(r, 'md')}\n\n"
            out = out.strip()

    if args.output_file:
        of = open(args.output_file, "w", encoding="utf-8")        
        of.write(out)
        of.close()
    else:
        print(out)

if __name__ == "__main__":
    main(sys.argv[1:])

# PubMed fields ###################################################################################
"""
## PubMed fields

### Most useful

[TI] Title  
[AB] Abstract  
[PT] Publication Type  
[DP] Date of Publication  
[LA] Language  
[PMID] PubMed Unique Identifier  
[MH] MeSH Terms  
[SO] Source  

#### Authors

Sequential sets for each author of  

[FAU] Full Author  
[AU] Author  
[AUID] Author Identifier (not always)  
[AD] Affiliation (multiple entries possible)  

#### Publication

[TA] Journal Title Abbreviation  
[JT] Journal Title  

[CI] Copyright Information  
[IRAD] Investigator Affiliation  
[AID] Article Identifier  
[BTI] Book Title  
[CTI] Collection Title  
[] Comments/Corrections	   
[COIS] Conflict of Interest Statement  
[CN] Corporate Author  
[CRDT] Create Date  
[DCOM] Date Completed  
[DA] Date Created  
[LR] Date Last Revised  
[DEP] Date of Electronic Publication  
[EN] Edition  
[ED] Editor  
[FED] Full Editor Name^(.*?)[ \t\n]+\((.*?)\)  
[EDAT] Entrez Date  
[GS] Gene Symbol  
[GN] General Note  
[GR] Grant Number  
[IR] Investigator Name  
[FIR] Full Investigator Name  
[ISBN] ISBN  
[IS] ISSN  
[IP] Issue  
[TA] Journal Title Abbreviation  
[JT] Journal Title  
[LID] Location Identifier  
[MID] Manuscript Identifier  
[MHDA] MeSH Date  
[JID] NLM Unique ID  
[RF] Number of References  
[OAB] Other Abstract  
[OCI] Other Copyright Information  
[OID] Other ID  
[OT] Other Term  
[OTO] Other Term Owner  
[OWN] Owner  
[PG] Pagination  
[PS] Personal Name as Subject  
[FPS] Full Personal Name as Subject  
[PL] Place of Publication  
[PHST] Publication History Status  
[PST] Publication Status  
[PUBM] Publishing Model  
[PMC] PubMed Central Identifier  
[PMCR] PubMed Central Release  
[RN] Registry Number/EC Number  
[NM] Substance Name  
[SI] Secondary Source ID  
[SFM] Space Flight Mission  
[STAT] Status  
[SB] Subset  
[TT] Transliterated Title  
[VI] Volume  
[VTI] Volume Title
"""
# PubMed format ###################################################################################
""" format
PMID- 31281835
OWN - NLM
STAT- MEDLINE
DCOM- 20191212
LR  - 20220409
IS  - 2314-6141 (Electronic)
IS  - 2314-6133 (Print)
VI  - 2019
DP  - 2019
TI  - Diagnostic Accuracy of Lever Sign Test in Acute, Chronic, and Postreconstructive 
      ACL Injuries.
PG  - 3639693
LID - 10.1155/2019/3639693 [doi]
LID - 3639693
AB  - BACKGROUND: The aim of this study is to determine the diagnostic accuracy of 
      lever sign test in acute, chronic, and postreconstructive ACL injuries. METHODS: 
      In total, 78 patients (69 male, 9 female) were subjected to clinical instability 
      tests including Lachman, anterior drawer, pivot shift, and lever sign when an 
      injury of the ACL was suspected. All tests were performed bilaterally in all 
      patients in acute, chronic period and patients who underwent surgery after the 
      anaesthesia and after the reconstruction at the last follow-up by two senior 
      orthopaedic surgeons. MRI was taken from all patients and MRI image was taken as 
      the reference test when evaluating the accuracy of the tests. RESULTS: The mean 
      age of patients was 26.2±6.4 years (range, 17-44 years). Sensitivity and accuracy 
      values of the Lachman, anterior drawer, pivot shift, and lever tests in the acute 
      phase were calculated as 80.6%, 77.4%, 51.6%, 91.9% and 76.9%, 75.6%, 60.3%, 
      92.3%, respectively, and in the chronic (preanaesthesia) phase were calculated as 
      83.9%, 79.0%, 56.5%, 91.9% and 80.8%, 78.2%, 64.1%, 92.3%, respectively. Lachman, 
      anterior drawer, pivot shift, and lever sign Acute's significant [AUC: 0.716, 
      0.731, 0.727, 0.928, respectively] activity were observed in the prediction of 
      ACL rupture in MRI. CONCLUSION: An ideal test to diagnose the integrity of the 
      ACL should be easy to perform and reproducible with high sensitivity and 
      specificity. From this perspective, the lever test seems to be a good test for 
      clinicians in acute, chronic and postreconstructive ACL injuries.
FAU - Gürpınar, Tahsin
AU  - Gürpınar T
AD  - Istanbul Training and Research Hospital, Department of Orthopedics and 
      Traumatology, Istanbul, Turkey.
FAU - Polat, Barış
AU  - Polat B
AUID- ORCID: 0000-0001-8229-6412
AD  - University of Kyrenia, Faculty of Medicine, Department of Orthopaedics and 
      Traumatology, Kyrenia, Cyprus.
FAU - Polat, Ayşe Esin
AU  - Polat AE
AD  - Dr. Akçiçek State Hospital, Department of Orthopaedics and Traumatology, Kyrenia, 
      Cyprus.
FAU - Çarkçı, Engin
AU  - Çarkçı E
AD  - Istanbul Training and Research Hospital, Department of Orthopedics and 
      Traumatology, Istanbul, Turkey.
FAU - Öztürkmen, Yusuf
AU  - Öztürkmen Y
AUID- ORCID: 0000-0002-2199-2411
AD  - Istanbul Training and Research Hospital, Department of Orthopedics and 
      Traumatology, Istanbul, Turkey.
LA  - eng
PT  - Journal Article
DEP - 20190609
PL  - United States
TA  - Biomed Res Int
JT  - BioMed research international
JID - 101600173
SB  - IM
MH  - Acute Disease
MH  - Adolescent
MH  - Adult
MH  - Anterior Cruciate Ligament/surgery
MH  - Anterior Cruciate Ligament Injuries/*diagnosis/*surgery
MH  - *Anterior Cruciate Ligament Reconstruction
MH  - Chronic Disease
MH  - *Diagnostic Tests, Routine
MH  - Female
MH  - Humans
MH  - Magnetic Resonance Imaging
MH  - Male
MH  - Meniscus/surgery
MH  - Physical Examination
MH  - Sensitivity and Specificity
MH  - Young Adult
PMC - PMC6590604
EDAT- 2019/07/10 06:00
MHDA- 2019/12/18 06:00
CRDT- 2019/07/09 06:00
PHST- 2019/01/24 00:00 [received]
PHST- 2019/04/28 00:00 [accepted]
PHST- 2019/07/09 06:00 [entrez]
PHST- 2019/07/10 06:00 [pubmed]
PHST- 2019/12/18 06:00 [medline]
AID - 10.1155/2019/3639693 [doi]
PST - epublish
SO  - Biomed Res Int. 2019 Jun 9;2019:3639693. doi: 10.1155/2019/3639693. eCollection 
      2019.
"""