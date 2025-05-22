using Documenter, TinyCrypto
push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "TinyCrypto.jl",
    modules = [TinyCrypto],
    pages = [
        "Home" => "index.md",
        "Utility Functions" => "utils.md",
    ],
)

deploydocs(
    repo = "github.com/vargaconsulting/tiny-crypto.git",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "",          # important: make 'dev' go to root
    target = "build",     # matches your makedocs target
    dirname = "."         # means push to root of gh-pages
)