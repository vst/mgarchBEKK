{ ... }:

let
  ## Import this codebase's Nix helper set:
  nix = import ./nix { };

  ## Get packages:
  pkgs = nix.pkgs;

  ## R environment:
  r = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [
      devtools
      mvtnorm
      roxygen2
      spelling
      testthat
      tseries
    ];
  };
in
pkgs.mkShell {
  buildInputs = [
    ## Fancy stuff:
    pkgs.figlet
    pkgs.lolcat

    ## Development:
    pkgs.texlive.combined.scheme-full
    r
  ];

  shellHook = ''
    figlet -w 999 "mgarchBEKK DEV SHELL" | lolcat -S 42
  '';
}
