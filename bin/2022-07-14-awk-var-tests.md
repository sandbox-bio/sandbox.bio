Currently, `awk -v abc=3 '{ print $2 }'` gives an error because shell-parse doesn't support variable assignment mid-command.

I tried modifying `grammar.pegjs` in `node_modules/shell-parse/` by removing this line:
```
!variableAssignment
```

Then I installed pegjs (`npm install pegjs`) and `npm install git://github.com/for-GET/pegjs-override-action.git#semver:~0.5`

Then I rebuilt `parser.js` using:

```bash
node ./build.js > node_modules/shell-parse/parser.js
```

The diff between the two:

```diff
index 48f4f91..b4c2535 100644
--- a/node_modules/shell-parse/parser.js.bak
+++ b/node_modules/shell-parse/parser.js
@@ -1,6 +1,6 @@
 module.exports = parse
 
-function parse(input, opts) {
+function parse (input, opts) {
   // Wrap parser.parse to allow specifying the start rule
   // as a shorthand option
   if (!opts) {
@@ -1964,7 +1964,7 @@ var parser=(function() {
     }
 
     function peg$parsecommandName() {
-      var s0, s1, s2, s3, s4;
+      var s0, s1, s2, s3;
 
       peg$silentFails++;
       s0 = peg$currPos;
@@ -1990,29 +1990,14 @@ var parser=(function() {
           s2 = peg$FAILED;
         }
         if (s2 !== peg$FAILED) {
-          s3 = peg$currPos;
-          peg$silentFails++;
-          s4 = peg$parsevariableAssignment();
-          peg$silentFails--;
-          if (s4 === peg$FAILED) {
-            s3 = void 0;
-          } else {
-            peg$currPos = s3;
-            s3 = peg$FAILED;
+          s3 = peg$parseconcatenation();
+          if (s3 === peg$FAILED) {
+            s3 = peg$parsebuiltinCommandName();
           }
           if (s3 !== peg$FAILED) {
-            s4 = peg$parseconcatenation();
-            if (s4 === peg$FAILED) {
-              s4 = peg$parsebuiltinCommandName();
-            }
-            if (s4 !== peg$FAILED) {
-              peg$savedPos = s0;
-              s1 = peg$c70(s4);
-              s0 = s1;
-            } else {
-              peg$currPos = s0;
-              s0 = peg$FAILED;
-            }
+            peg$savedPos = s0;
+            s1 = peg$c70(s3);
+            s0 = s1;
           } else {
             peg$currPos = s0;
             s0 = peg$FAILED;
@@ -3935,22 +3920,22 @@ var parser=(function() {
 
     var isArray = require("isarray")
     var map = require("array-map")
-    function join(arr) {
+    function join (arr) {
         return arr.join("")
       }
-    function literal(string) {
+    function literal (string) {
         return {
           type: 'literal',
           value: isArray(string) ? string.join('') : string
         }
       }
-    function first(arr) {
+    function first (arr) {
         return arr[0]
       }
-    function second(arr) {
+    function second (arr) {
         return arr[1]
       }
-    function flattenConcatenation(pieces) {
+    function flattenConcatenation (pieces) {
         // TODO - this algo sucks, it's probably on the order of n^4
         var result = [pieces[0]]
           , len = pieces.length
```

Although that made `awk -v abc=3 '{ print $2 }'` work, now variable assignments like `ABC=123` stopped working :(


