The used functions do not behave equal for all parameter settings.
Examples are provided for the commutativity and the transitive law. This
affects for example the cardinality of the joined sets. Moreover, if the
distribution of the combinations within a set is studied, it should be
taken into account if each feature is counted or each combination
(e.g. *A* = {(*a*<sub>1</sub>, *b*<sub>1</sub>), (*a*<sub>1</sub>, *b*<sub>2</sub>)}).

## Laws to consider for bedtools closest

    cat example/A.bed

    chr1    50  60
    chr1    100 200

    cat example/B.bed

    chr1    100 200
    chr1    200 300

### Commutative property

The nearest neighbor function is not commutative for all parameter
settings.

    bedtools closest -D a -a example/A.bed -b example/B.bed

    chr1    50  60  chr1    100 200 41
    chr1    100 200 chr1    100 200 0

    bedtools closest -D a -a example/B.bed -b example/A.bed

    chr1    100 200 chr1    100 200 0
    chr1    200 300 chr1    100 200 -1

Therefore it is necessary to consider the symmetry of the used
parameters. For example by using `ignore overlap`, `ignore downstream`
and `ignore upstream`,

    bedtools closest -D a -io -iu -a example/A.bed -b example/B.bed

    chr1    50  60  chr1    100 200 41
    chr1    100 200 chr1    200 300 1

    bedtools closest -D a -io -id -a example/B.bed -b example/A.bed

    chr1    100 200 chr1    50  60  -41
    chr1    200 300 chr1    100 200 -1

or by enumerating the complete environment of a feature.

    bedtools closest -D a -k 2 -a example/A.bed -b example/B.bed

    chr1    50  60  chr1    100 200 41
    chr1    50  60  chr1    200 300 141
    chr1    100 200 chr1    100 200 0
    chr1    100 200 chr1    200 300 1

    bedtools closest -D a -k 2 -a example/B.bed -b example/A.bed

    chr1    100 200 chr1    100 200 0
    chr1    100 200 chr1    50  60  -41
    chr1    200 300 chr1    100 200 -1
    chr1    200 300 chr1    50  60  -141
