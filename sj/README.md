# Stata Journal submission materials

`phdid_sj.tex` is the draft article.

## Typesetting

The Stata Journal class files (`sj.cls`, `stata.sty`, and the `stlog`
environment support) are distributed by the Journal and are deliberately not
vendored here. Download them from

  https://www.stata-journal.com/submissions/

place them beside `phdid_sj.tex`, and run `pdflatex phdid_sj` twice.

## What the Journal asks for, and where it is

| Requirement | Where |
|---|---|
| Article source | `sj/phdid_sj.tex` |
| Ado-files | `ado/` |
| Help files | `doc/*.sthlp` |
| Certification script | `cert/phdid_cert.do` |
| Worked example that reproduces the article's output | `examples/mpdta.do` |
| Package files for `net install` | `phdid.pkg`, `stata.toc` |

## Before submitting

1. `do cert/phdid_cert.do` — must report all checks passed.
2. `do examples/mpdta.do` — must run clean and produce the numbers quoted in
   section 4 of the article.
3. Regenerate any figures included in the article from `examples/mpdta.do`.
