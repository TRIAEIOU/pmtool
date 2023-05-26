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