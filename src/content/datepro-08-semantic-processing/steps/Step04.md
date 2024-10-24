<script>
import Execute from "$components/Execute.svelte";
</script>

##  Ancestors

Finding all the ancestors of a class includes many chain invocations of the `getparents.sh` until we get no matches. We also should avoid relations that
are cyclic, otherwise we will enter in a infinite loop. Thus, for identifying the ancestors of a class, we will only consider parent relations, i.e. subsumption relations.

#### Grandparents

In the previous step we were able to extract the direct parents of a class, but the parents of these parents also represent generalizations of the original class. For example, to get the parents of the parents (grandparents) of
malignant hyperthermia we need to invoke getparents.sh twice:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl | ./getparents.sh doid.owl | ./getparents.sh doid.owl" />

And we will find the URIs of the grandparents of malignant hyperthermia.

Or to get their labels we can add the `getlabels.sh` script:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl | ./getparents.sh doid.owl | ./getparents.sh doid.owl | ./getlabels.sh doid.owl" />

And we find the labels of the grandparents of malignant hyperthermia.

#### Root class
However, there are classes that do not have any parent, which are called root classes. _disease_ and _chemical entity_ are root classes of DO and ChEBI ontologies, respectively. As we can see these are
highly generic terms.
To check if it is the root class, we can ask for their parents:

<Execute command="echo 'disease' | ./geturi.sh doid.owl | ./getparents.sh doid.owl" />

<Execute command="echo 'chemical entity' | ./geturi.sh chebi_lite.owl | ./getparents.sh chebi_lite.owl" />

In both cases, we will get the warning that no matches were found, confirming that they are the root class.
```text
XPath set is empty
```

#### Recursion
We can now build a script that receives a list of URIs as standard input, and invokes `getparents.sh` recursively until it reaches the root class. The script named `getancestors.sh`

<Execute command="nano getancestors.sh" />

should contain the following lines:

```bash
OWLFILE=$1
CLASSES=$(cat -)
[[ -z "$CLASSES" ]] && exit
PARENTS=$(echo "$CLASSES" | ./getparents.sh $OWLFILE | sort -u)
echo "$PARENTS"
echo "$PARENTS" | ./getancestors.sh $OWLFILE
```

The second line of the script saves the standard input in a variable named `CLASSES`, because we need to use it twice: i) to check if the input as any classes or is empty (line 3) and ii) to get the parents of the classes given as input (line 4). If the input is empty then the script ends, this is the base case
of the [recursion](https://en.wikipedia.org/wiki/Recursion). This is required so the recursion stops at a given point. Otherwise, the script would run indefinitely until the user stops it manually. The fourth line of the script stores the output in a variable named `PARENTS`, because we need also to use it twice: i) to output these direct parents (line 5), and ii) to get the ancestors of these parents (line 6). We should note that we are invoking the `getancestors.sh` script inside the `getancestors.sh`, which
defines the recursion step. Since the subsumption relation is acyclic, we expect that at some time we will reach classes without parents (root classes) and then the script will end.
We should note that the echo of the variables `CLASSES` and `PARENTS` need to be inside commas, so the newline characters are preserved.


#### Iteration
Recursion is most of the times computational expensive, but usually it is possible to replace recursion with iteration to develop a more efficient algorithm.
Explaining iteration and how to refactor a recursive script is out of scope of this tutorial, nevertheless the following script represents an equivalent way to get all the ancestors without using recursion:

```bash
OWLFILE=$1
CLASSES=$(cat -)
ANCESTORS=""
while [[ ! -z "$CLASSES" ]]
do
  PARENTS=$(echo "$CLASSES" | ./getparents.sh $OWLFILE | sort -u)
  ANCESTORS="$ANCESTORS\n$PARENTS"
  CLASSES=$PARENTS
done
echo -e "$ANCESTORS"
```

The script uses the while command that basically implements iteration by repeating a set of commands (lines 6-8) while a given condition is satisfied (line 4).

To test the recursive script, we can provide as standard input the label malignant hyperthermia:

<Execute command="chmod u+x getancestors.sh" />

<Execute command="echo 'http://purl.obolibrary.org/obo/DOID_8545' | ./getancestors.sh doid.owl" />

The output will be the URI of all its ancestors

We should note that we will still receive the XPath warning when the script reaches the root class and no parents are found:
```text
XPath set is empty
```

To remove this warning and just get the labels of the ancestors of malignant hyperthermia, we can redirect the warnings to the null device:

<Execute command="echo 'malignant hyperthermia' | ./geturi.sh doid.owl | ./getancestors.sh doid.owl 2>/dev/null | ./getlabels.sh doid.owl" />

The output will now include the name of all ancestors of malignant hyperthermia.
We should note that the first two ancestors are the direct parents of malignant hyperthermia, and the last one is the root class. This happens because the recursive script prints the parents before invoking itself to find the ancestors
of the direct parents.
We can do the same with caffeine, but be advised that given the higher number of ancestors in ChEBI we may now have to wait a little longer for the script to end.

<Execute command="echo 'caffeine' | ./geturi.sh chebi_lite.owl | ./getancestors.sh chebi_lite.owl | ./getlabels.sh chebi_lite.owl | sort -u" /> 

The results include repeated classes that were found by using different branches, so that is why we need to add the sort command with the `-u` option to eliminate the duplicates. The script will print the ancestors being found by the script.