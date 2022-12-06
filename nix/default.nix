{ ...
}:

let
  ## Import sources:
  sources = import ./sources.nix;

  ## Import nixpkgs:
  pkgs = import sources.nixpkgs { };

  ## Define R package dependencies:
  rDeps = with pkgs.rPackages; {
    ## Dependencies for our package:
    pkg = [
      mvtnorm
      tseries
    ];

    ## Dependencies for development:
    dev = [
      autonewsmd
      devtools
      roxygen2
      spelling
      testthat
    ];
  };

  ## Setup R environment:
  r = pkgs.rWrapper.override {
    packages = rDeps.pkg ++ rDeps.dev;
  };

  ## Get aspell with dictionaries:
  aspell = pkgs.aspellWithDicts (d: [ d.en ]);

  ## Setup development shell:
  shell = pkgs.mkShell
    {
      buildInputs = [
        ## Fancy stuff:
        pkgs.figlet
        pkgs.lolcat

        ## Spell checking related:
        aspell

        ## TeX stuff:
        pkgs.texlive.combined.scheme-full

        ## R:
        r
      ];

      shellHook = ''
        figlet -w 999 "mgarchBEKK DEV SHELL" | lolcat -S 42

        export ASPELL_CONF="dict-dir ${aspell}/lib/aspell"
      '';
    };
in
{
  pkgs = pkgs;
  r = r;
  shell = shell;
}
