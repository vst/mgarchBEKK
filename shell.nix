{ ... }:

let
  ## Import this codebase's Nix helper set:
  nix = import ./nix { };
in
nix.shell
