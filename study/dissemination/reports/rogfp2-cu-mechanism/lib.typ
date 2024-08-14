//////////////////////
// Helper functions //
//////////////////////

#let todo(some_text, color: rgb("ff8c00"), weight: "medium") = {
    text(fill: color, weight: "black", "TODO: ")
    text(fill: color, weight: weight, some_text)
}

#let addcite(info: "") = {
    let str = "(CITE"
    if info != "" {
        str += " " + info
    }
    str += ")"
    text(fill: rgb("5203fc"), weight: "black", str)
}

////////////////
// Formatting //
////////////////

#let std-bibliography = bibliography

#let config(
    // The paper's title.
    title: [Paper Title],

    // An array of authors. For each author you can specify a name,
    // department, organization, location, and email. Everything but
    // but the name is optional.
    authors: (),
    version: "",
    bibliography: none,
    body
) = {

    show link: set text(fill: rgb("#0077b6"))

    // Set document metadata.
    set document(title: title, author: authors.map(author => author.name))

    // Set the body font.
    set text(font: "Roboto", lang: "en", size: 10pt, hyphenate: false)

    // Tables & figures
    set figure(placement: top)
    show figure.where(kind: table): set figure.caption(position: top)
    show figure.where(kind: table): set text(size: 8pt)
    show figure.caption.where(kind: table): smallcaps
    show figure.where(kind: table): set figure(numbering: "I")

    show figure.where(kind: image): set figure(supplement: [Fig.], numbering: "1")
    show figure.caption: set text(size: 8pt)

    // Code blocks
    show raw: set text(font: "TeX Gyre Cursor", size: 10pt)

    // Configure the page.
    set page(
        paper: "us-letter",
        margin: (left: 0.75in, right: 0.75in, top: 0.75in, bottom: 0.75in),
            footer: [
            #set text(9pt)

            #grid(
                columns: (50%, 50%),
                rows: (1em),
                align(
                    left,
                    text(
                        weight: "regular"
                    )[#version]
                ),
                align(
                    right,
                    counter(page).display(
                        "1 of 1",
                        both: true,
                    )
                ),
            )
        ],
    )

    // Configure headings
    set heading(numbering: "1.A.1.")
    show heading: it => locate(loc => {
        // Find out the final number of the heading counter.
        let levels = counter(heading).at(loc)
        let deepest = if levels != () {
            levels.last()
        } else {
            1
        }

        set text(12pt, weight: 400)
        if it.level == 1 [
            #set text(size: 12pt, weight: "medium")
            #set align(center)
            #v(20pt, weak: true)
            #if it.numbering != none{
                numbering("1.", deepest)
                h(4pt, weak: true)
            }
            #it.body
            #v(13.75pt, weak: true)
        ] else if it.level == 2 [
            // Second-level headings are run-ins.
            #set par(first-line-indent: 0pt)
            #set text(style: "italic")
            #v(10pt, weak: true)
            #if it.numbering != none {
                numbering("A.", deepest)
                h(7pt, weak: true)
            }
            #it.body
            #v(10pt, weak: true)
        ] else [
            // Third level headings are run-ins too, but different.
            #if it.level == 3 {
                numbering("a.", deepest)
                [ ]
            }
            _#(it.body):_
        ]
    })

    // Display the paper's title.
    v(3pt, weak: true)
    align(center, text(18pt, title))
    v(8.35mm, weak: true)

    // Display the authors list.
    for i in range(calc.ceil(authors.len() / 3)) {
        let end = calc.min((i + 1) * 3, authors.len())
        let is-last = authors.len() == end
        let slice = authors.slice(i * 3, end)
        grid(
        columns: slice.len() * (1fr,),
        gutter: 12pt,
        ..slice.map(author => align(center, {
            text(12pt, author.name)
            if "department" in author [
            \ #emph(author.department)
            ]
            if "organization" in author [
            \ #emph(author.organization)
            ]
            if "location" in author [
            \ #author.location
            ]
            if "email" in author [
            \ #link("mailto:" + author.email)
            ]
        }))
        )

        if not is-last {
            v(16pt, weak: true)
        }
    }
    v(40pt, weak: true)

    // Start two-column mode
    show: columns.with(2, gutter: 12pt)

    outline(indent: 1em, depth: 1)

    set par(justify: true, first-line-indent: 1em)
    show par: set block(spacing: 0.65em)

    // Display the paper's contents.
    body

    // Display bibliography.
    show heading: it => [
        #set par(justify: true, first-line-indent: 0em)
        #v(0.5em)
        #it.body
        #v(0.5em)
    ]
    if bibliography != none {
        show std-bibliography: set text(10pt)
        bibliography
    }
}
