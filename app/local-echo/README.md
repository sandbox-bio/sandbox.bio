# local-echo

Copied from [wavesoft/local-echo](https://github.com/wavesoft/local-echo) repo; commit `8d0b7f5`

### Changes

Changes are denoted by `FIXED:` labels

* History: Limit removes the newest command, not the oldest; [7e5e3e1](https://github.com/sandbox-bio/sandbox.bio/commit/7e5e3e1)
* History: Fix bug where if ran the previous command again, the history pointer would shift, which is very confusing; [3b9f5b0](https://github.com/sandbox-bio/sandbox.bio/commit/3b9f5b0)
* Cursor: When entering multiline inputs, the cursor is off by the number of breaklines; [f4cd9bd](https://github.com/sandbox-bio/sandbox.bio/commit/f4cd9bd)
* Cursor: Make sure cursor behaves well with a color prompt	; [3c3155a](https://github.com/sandbox-bio/sandbox.bio/commit/3c3155a)
