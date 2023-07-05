<p class="lead">sandbox.bio features interactive tutorials for exploring bioinformatics command-line tools in a secure sandbox.</p>

<p class="lead fw-bold mt-4 mb-1">How does it work?</p>
<p>We compiled commonly-used bioinformatics tools (and general tools such as <code>awk</code>, <code>grep</code> and <code>coreutils</code>!) to WebAssembly so that they can run in your browser (see our <a href="https://github.com/biowasm/biowasm" target="_blank" rel="noreferrer">biowasm</a> project).</p>

<p>To simulate a terminal environment, command-line input is parsed into an abstract syntax tree and supports features such as piping (<code>|</code>), file redirection (<code>></code>, <code>>></code>), process substitution (<code>&lt;()</code>), conditional commands (<code>&&</code>, <code>||</code>), asynchronous commands (<code>&</code>), variables (<code>abc=123</code>), and autocomplete. The terminal UI is based on <a href="https://github.com/xtermjs/xterm.js/" target="_blank">xterm.js</a>. All the files you read and write to are temporarily stored in memory using Emscripten's <a href="https://emscripten.org/docs/api_reference/Filesystem-API.html" target="_blank" rel="noreferrer">virtual file system</a>.</p>

<p>The <a href="/playgrounds?id=awk">playgrounds</a> use sample data from <a href="https://github.com/TheUpshot/chipotle">The UpShot</a>.</p>

<p class="lead fw-bold mt-4 mb-1">How to contribute</p>
<p>If you have feedback, ideas for new tutorials, or if one of your own bioinformatics tutorials could benefit from being interactive, please <a href="https://github.com/sandbox-bio/sandbox.bio/discussions" target="_blank" rel="noreferrer">reach out</a>!</p>

<p class="lead fw-bold mt-4 mb-1">Author</p>
<p>Built by <a href="https://robert.bio" target="_blank" rel="noreferrer">Robert Aboukhalil</a>.</p>
