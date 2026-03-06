<img src="docs/img/peppro_logo.svg" alt="pepatac logo" height="200" align="left"/>  

<br></br>
<br></br>
<br></br>
---

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pep.databio.org)

PEPPRO is a pipeline designed to process PRO-seq (and GRO-seq) data. For more information see: http://peppro.databio.org/

## Install

```bash
pip install piper pipestat looper
```

**Note:** The pypiper PyPI package is `piper` (not `pypiper`, which is an unrelated package).

## Testing

Unit tests need no special setup:

```bash
pytest tests/test_unit.py -v
```

Integration tests require bioinformatics tools via [bulker](https://bulker.io). Use the wrapper script:

```bash
bash tests/scripts/test-integration.sh
```

This runs `bulker exec databio/peppro:1.1.0` to provide samtools, bowtie2, bedtools, etc. Do NOT run integration tests without bulker — they will fail with missing tools.

## Docs

Develop docs with:

```
mkdocs serve -f mkdocs.yml
```

Build and deploy with:

```
mkdocs build -f mkdocs.yml -d $CODEBASE/code.databio.org/peppro
```

## Contributing

Pull requests welcome. Active development should occur in a development or feature branch.

## Contributors

* Jason Smith, jasonsmith@virginia.edu
* Michael Guertin, mjg7y@virginia.edu
* Nathan Sheffield, nathan@code.databio.org
* Others... (add your name)

