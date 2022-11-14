{ ... }:

let
  ## Import this codebase's Nix helper set:
  nix = import ./nix { };

  ## Get packages:
  pkgs = nix.pkgs;

  ## R environment:
  r = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [
      tseries
      mvtnorm
      testthat
      devtools
      roxygen2
    ];
  };
in
pkgs.mkShell {
  buildInputs = [
    ## Fancy stuff:
    pkgs.figlet
    pkgs.lolcat

    ## Development:
    r
  ];

  shellHook = ''
    figlet -w 999 "mgarchBEKK DEV SHELL" | lolcat -S 42
  '';
}
