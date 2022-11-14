{ compiler ? "ghc90"
, ...
}:

let
  ## Import sources:
  sources = import ./sources.nix;

  ## Import telosnix:
  telosnix = import sources.telosnix { };

  ## Import nixpkgs:
  pkgs = import telosnix.pkgs-sources.unstable { };
in
{
  sources = sources;
  telosnix = telosnix;
  pkgs = pkgs;
}
