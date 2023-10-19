# Standard operating procedures

## Markdown files

All markdown files should follow the following formatting procedures:

- One complete sentence per line, regardless of length;
- Indents shall be two spaces;
- Use [MyST-flavored](https://myst-parser.readthedocs.io/en/latest/) formatting.

Each Markdown file should have a `toctree` at the bottom that allows our website framework to find it.
For example, the `index.md` page has

````markdown
```{toctree}
:hidden:

website/source/environment
01-management/README
02-methods/README
03-data/README
04-analysis/README
05-figures/README
06-dissemination/README
07-misc/README
website/source/sop
Code license <LICENSE_CODE>
Info license <LICENSE_INFO>
```
````

## Obsidian

It is **highly recommended** to make the repository root an [Obsidian](https://obsidian.md) vault.

We recommend the following plugins:

- [obsidian-calendar-plugin](https://github.com/liamcain/obsidian-calendar-plugin),
- [cMenu-plugin](https://github.com/chetachiezikeuzor/cMenu-Plugin),
- [obsidian-dataview](https://github.com/blacksmithgu/obsidian-dataview),
- [obsidian-table-editor](https://github.com/ganesshkumar/obsidian-table-editor),
- [obsidian-bibtex-adder](https://github.com/oasci/obsidian-bibtex-adder),
- [obsidian-citation-plugin](https://github.com/hans/obsidian-citation-plugin),
- [obsidian-frontmatter-tag-suggest](https://github.com/jmilldotdev/obsidian-frontmatter-tag-suggest),
- [obsidian-git](https://github.com/denolehov/obsidian-git),
- [obsidian-open-link-with](https://github.com/MamoruDS/obsidian-open-link-with),
- [obsidian-paste-mode](https://github.com/jglev/obsidian-paste-mode),
- [obsidian-table-editor](https://github.com/ganesshkumar/obsidian-table-editor),
- [tag-wrangler](https://github.com/pjeby/tag-wrangler).

## Website

We use [sphinx](https://www.sphinx-doc.org/en/master/) to automatically generate our website using the [furo](https://github.com/pradyunsg/furo) theme.
For guidelines of what formatting options you have, please see [MyST](https://myst-parser.readthedocs.io/en/latest/).

To locally view the website (assuming you have [setup the environment](./environment)) run the following command.

```bash
make website
```
