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

  ## Setup development shell:
  shell = pkgs.mkShell
    {
      buildInputs = [
        ## Fancy stuff:
        pkgs.figlet
        pkgs.lolcat

        ## TeX stuff:
        pkgs.texlive.combined.scheme-full

        ## R:
        r
      ];

      shellHook = ''
        figlet -w 999 "mgarchBEKK DEV SHELL" | lolcat -S 42
      '';
    };
in
{
  pkgs = pkgs;
  r = r;
  shell = shell;
}
