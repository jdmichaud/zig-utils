# Zig utils

Requires Zig 0.17. `clap.zig` and `geometry.zig` build with it; the other modules
still need porting from Zig 0.14.

- `clap.zig`: Command line argument parser
- `draw.zig`: Provide a canvas context 'like' interface for drawing
- `io-adapter`:
  - Keyboard/mouse input + graphic output
  - One adapter for SDL for now
- `misc.zig`:
  - Load a file
  - Write an image (pgm or pam)

