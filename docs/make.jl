using Documenter

import PALEOcopse

using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src/paleo_references.bib"))

makedocs(;
    sitename = "PALEOcopse Documentation", 
    pages = [
        "index.md",
        "Examples and Tutorials" => [
            "ExampleCOPSEreloaded.md",
        ],
        "Design" => [
            "COPSE_Domains.md",
        ],
        # no HOWTO docs yet
        "Reference" => [
            "PALEOcopse.md",
        ],
        "References.md",
        "indexpage.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    plugins = [bib],
)

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOcopse.jl.git",
)
