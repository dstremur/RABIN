{
  description = "RABIN dev shell";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-26.05";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in
    {
      devShells.${system}.default = pkgs.mkShell {
        packages = with pkgs; [
          gmp
          gcc
		  gnumake 
          flint
          doxygen
          openssl
        ];

        shellHook = ''
          export NIX_ENFORCE_NO_NATIVE=0
        '';
      };
    };
}
