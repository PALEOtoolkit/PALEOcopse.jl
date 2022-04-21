using Documenter

import PALEOcopse

using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src/paleo_references.bib"))

makedocs(bib, sitename="PALEOcopse Documentation", 
# makedocs(sitename="PALEO Documentation", 
    pages = [
        "index.md",
        "Examples and Tutorials" => [
            "ExampleCOPSEreloaded.md",
        ],
        # no Design docs yet
        # no HOWTO docs yes
        "Reference" => [
            "PALEOcopse.md",
        ],
        "References.md",
        "indexpage.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOcopse.jl.git",
)
