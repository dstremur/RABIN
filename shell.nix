{ pkgs ? import <nixpkgs> {}}:

pkgs.mkShell {
  buildInputs = [
    pkgs.gmp
    pkgs.gcc
    pkgs.flint
    pkgs.doxygen
    pkgs.openssl
  ];
}
