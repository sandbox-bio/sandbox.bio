# local-echo

Copied from [wavesoft/local-echo](https://github.com/wavesoft/local-echo) repo; commit `8d0b7f5`

### Changes

Changes are denoted by `FIXED:` labels

- History: Limit removes the newest command, not the oldest; [7e5e3e1](https://github.com/sandbox-bio/sandbox.bio/commit/7e5e3e1)
- History: Fix bug where if ran the previous command again, the history pointer would shift, which is very confusing; [3b9f5b0](https://github.com/sandbox-bio/sandbox.bio/commit/3b9f5b0)
- ~~Cursor: When entering multiline inputs, the cursor is off by the number of breaklines; [f4cd9bd](https://github.com/sandbox-bio/sandbox.bio/commit/f4cd9bd)~~ Undone by commit [221eb39](https://github.com/sandbox-bio/sandbox.bio/commit/221eb39)
- Cursor: Make sure cursor behaves well with a color prompt; [3c3155a](https://github.com/sandbox-bio/sandbox.bio/commit/3c3155a)
- Cursor: When entering multiline inputs, the cursor is off (attempt 2); [1c3dcd5](https://github.com/sandbox-bio/sandbox.bio/commit/1c3dcd5)
- Cursor: Using "\" with a color prompt creates extra whitespace; [2f2248e](https://github.com/sandbox-bio/sandbox.bio/commit/2f2248e)

### ANSI Codes

See https://invisible-island.net/xterm/ctlseqs/ctlseqs.html

```
ESC A     Cursor up
ESC B     Cursor down
ESC C     Cursor right
ESC D     Cursor left
```
