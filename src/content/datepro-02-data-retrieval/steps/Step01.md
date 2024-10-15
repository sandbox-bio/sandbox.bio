<script>
import Execute from "$components/Execute.svelte";
</script>

The input argument(s) of our retrieval task is the chemical compound(s) of
which we want to retrieve more information. For the sake of simplicity, we
will start by assuming that the user knows the ChEBI identifier(s), i.e. the
script does not have to search by the name of the compounds.
So, the first step, is to automatically retrieve all proteins associated to the
given input chemical compound, that in our example is caffeine [CHEBI:27732](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27732).

#### Web Identifiers

In the manual process described in the book (Section: Caffeine Example), we downloaded the files by manually clicking on the links shown as Export options, namely the URLs:

1. https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=1&6578706f7274=1&chebiId=27732&dbName=UniProt
2. https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=2&6578706f7274=1&chebiId=27732&dbName=UniProt
3. https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-e=3&6578706f7274=1&chebiId=27732&dbName=UniProt

for downloading a CSV, Excel, or XML file, respectively.

We should note that the only difference between the three URLs is a single
numerical digit (1,` 2, and 3) after the first equals character (=), which means
that this digit can be used as an argument to select the type of file. Another
parameter that is easily observable is the ChEBI identifier (27732). Try to
replace 27732 by 17245 in any of those URLs, for
example:

- https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d-1169080-=1&6578706f7274=1&chebiId=17245&dbName=UniProt

Now we can use this new URL in the internet browser, and check what happens. If we did it correctly, our browser downloaded a file with more than
seven hundred proteins, since the 17245 is the ChEBI identifier of a popular
chemical compound in life systems, the carbon monoxide.

In this case, we are not using a fully RESTful web service, but the data path is pretty modular and self-explanatory. The path is clearly composed of:

- the name of the database (chebi);
- the method (viewDbAutoXrefs.do);
- and a list of parameters and their value (arguments) after the question
  mark character (?).

The order of the parameters in the URL is normally not relevant. They are
separated by the ampersand character (&) and the equals character (=) is
used to assign a value to each parameter (argument). This modular structure
of these URLs allows us to use them as data pipelines to fill our local files with
data, like pipelines that transport oil or gas from one container to another.

#### Single and double quotes

To construct the URL for a given ChEBI identifier, let us first understand the
difference between single quotes and double quotes in a string (sequence of
characters). We can create a script file named `getproteins.sh` by using a text
editor:

<Execute command="nano getproteins.sh" />

to add the following lines:

```bash
echo 'The input: $1'
echo "The input: $1"
```

Exit from the editor with **Ctrl-X**, and then press **Y** and **Enter** to save the file.

The command line tool `echo` displays the string received as argument.

Do not forget to save it in our working directory:

<Execute command="cat getproteins.sh" />

and add the right permissions with `chmod` as we did previously with our first script:

<Execute command="chmod u+x getproteins.sh" />

Now to execute the script we will only need to type:

<Execute command="./getproteins.sh" />

The output on the terminal should be:

```text
The input: $1
The input:
```

This means that when using single quotes, the string is interpreted literally
as it is, whereas the string within double quotes is analyzed, and if there is a
special character, such as the dollar sign (`$`), the script translates it to what
it represents. In this case, `$1` represents the first input argument. Since no
argument was given, the double quotes displays nothing.

To execute the script with an argument, we can type:

<Execute command="./getproteins.sh 27732" />

The output on our terminal should be:

```text
The input: $1
The input: 27732
```

We can check now that when using double quotes `$1` is translated to the
string given as argument.
Now we can update our script file named getproteins.sh

<Execute command="nano getproteins.sh" />

to contain only the following line:

```bash
echo "https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d
-1169080-e=1&6578706f7274=1&chebiId=$1&dbName=
UniProt"
```

Exit from the editor with **Ctrl-X**, and then press **Y** and **Enter** to save the file.

#### Comments

Instead of removing the previous lines, we can transform them in comments
by adding the hash character (`#`) to the beginning of the line:

```bash
#echo 'The input: $1'
#echo "The input: $1"
echo "https://www.ebi.ac.uk/chebi/viewDbAutoXrefs.do?d
-1169080-e=1&6578706f7274=1&chebiId=$1&dbName=
UniProt"
```

Commented lines are ignored by the computer during script execution, so this script will perform the same actions as the previous one with just one line of code.

Now, we can execute the script giving the ChEBI identifier as argument:

<Execute command="./getproteins.sh 2773" />

The output on our terminal should be the link that returns the CSV file containing the proteins associated with caffeine.

The next step will use this link to retrieve the data.
