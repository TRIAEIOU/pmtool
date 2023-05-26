# pmtool

Small python command line tool to query and/or parse PubMed from/to stdin or file. Output options are markdown (`md`) or (JSON `json`)

## Use

Command line options:

- `-i`/`--input-file`: read input from specified file (rather than stdin)
- `-o`/`--output-file`: write output to file (rather than stdout)
- `-f`/`--format`: output format (md/json, default md). If not specified but output file is, guess from filename
- `-n`/`--number`: max number of entries, `-1` for no limit (mainly useful for queries) (default -1)
- `-q`/`--query`: interpret input as queries to run against PubMed, remainder of command line interpreted as single query

Examples:

- Parse saved PubMed results (saved in `PubMed` format) to markdown file: `python pmtool -i saved-result-file-in-pubmed-format.txt -o formated-file.md`
- Parse saved PubMed results to stdout in json: `python pmtool -i saved-result-file-in-pubmed-format.txt -f json`
- Query PubMed: `python pmtool -o result.json -q some[ti] query[ab]`
- Run multiple queries retrieving max 100 results on each: `(echo query1 & echo Newline & echo query2) | python pmtool -o result.json -n 100 -q`
