{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {

  buildInputs = with pkgs; [
    gmp
	flint
	doxygen
	openssl
  ];

}
