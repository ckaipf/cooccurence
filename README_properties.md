The used functions do not behave equal for all parameter settings.
Examples are provided for the commutativity and the transitive law. This
affects for example the cardinality of the joined sets. Moreover, if the
distribution of the combinations within a set is studied, it should be
taken into account if each feature is counted or each combination
(e.g. *A* = {(*a*<sub>1</sub>, *b*<sub>1</sub>), (*a*<sub>1</sub>, *b*<sub>2</sub>)}).

## Laws to consider for bedtools closest

    cat A.bed

    chr1    50  60
    chr1    100 200

    cat B.bed

    chr1    100 200
    chr1    200 300

### Commutative property

The nearest neighbor function is not commutative for all parameter
settings.

    bedtools closest -D a -a A.bed -b B.bed

    chr1    50  60  chr1    100 200 41
    chr1    100 200 chr1    100 200 0

    bedtools closest -D a -a B.bed -b A.bed

    chr1    100 200 chr1    100 200 0
    chr1    200 300 chr1    100 200 -1

Therefore it is necessary to consider the symmetry of the used
parameters. For example by using `ignore overlap`, `ignore downstream`
and `ignore upstream`,

    bedtools closest -D a -io -iu -a A.bed -b B.bed

    chr1    50  60  chr1    100 200 41
    chr1    100 200 chr1    200 300 1

    bedtools closest -D a -io -id -a B.bed -b A.bed

    chr1    100 200 chr1    50  60  -41
    chr1    200 300 chr1    100 200 -1

or by enumerating the complete environment of a feature.

    bedtools closest -D a -k 2 -a A.bed -b B.bed

    chr1    50  60  chr1    100 200 41
    chr1    50  60  chr1    200 300 141
    chr1    100 200 chr1    100 200 0
    chr1    100 200 chr1    200 300 1

    bedtools closest -D a -k 2 -a B.bed -b A.bed

    chr1    100 200 chr1    100 200 0
    chr1    100 200 chr1    50  60  -41
    chr1    200 300 chr1    100 200 -1
    chr1    200 300 chr1    50  60  -141

## Properties and behaviour of dplyr::full\_join

### Commutative property

A non-commutative parameter setting may lead to inconsistent results. In
the second example, we used a not suitable parameter setting. The
obtained set *A* is then ambiguous: *a*<sub>2</sub> → *b*<sub>2</sub>
and *a*<sub>2</sub> → *N**A*.

    AB <- data.frame(A = c("a1", "a2"), B = c("b1", "b2"))
    BA <- data.frame(B = c("b1", "b2"), A = c("a1", "a2"))
    full_join(AB, BA)

    Joining, by = c("A", "B")

       A  B
    1 a1 b1
    2 a2 b2

    AB <- data.frame(A = c("a1", "a2"), B = c("b1", "b2"))
    BA <- data.frame(B = c("b1", NA), A = c("a1", "a2"))
    full_join(AB, BA)

    Joining, by = c("A", "B")

       A    B
    1 a1   b1
    2 a2   b2
    3 a2 <NA>

## Transitive relation

The same applies to transitivity: if (*a*<sub>1</sub>, *b*<sub>1</sub>)
and (*b*<sub>1</sub>, *c*<sub>1</sub>) were determined as joint
occurred, (*a*<sub>1</sub>, *c*<sub>1</sub>) should be implied. If that
is not true, multiple combinations are obtained in the union.

    AB <- data.frame(A = c("a1", "a2"), B = c("b1", "b2"))
    BA <- data.frame(B = c("b1", "b2"), A = c("a1", "a2"))
    AC <- data.frame(A = c("a1", "a2"), C = c("c1", "c2"))
    CA <- data.frame(C = c("c1", "c2"), A = c("a1", "a2"))

    BC <- data.frame(B = c("b1", "b2"), C = c("c1", "c2"))
    CB <- data.frame(C = c("c1", "c2"), B = c("b1", "b2"))

    AB_ <- full_join(AB, BA)

    Joining, by = c("A", "B")

    AC_ <- full_join(AC, CA)

    Joining, by = c("A", "C")

    A <- full_join(AB_, AC_)

    Joining, by = "A"

    BC_ <- full_join(BC, CB)

    Joining, by = c("B", "C")

    BA_ <- full_join(BA, AB)

    Joining, by = c("B", "A")

    B <- full_join(BA_, BC_)

    Joining, by = "B"

    full_join(A, B)

    Joining, by = c("A", "B", "C")

       A  B  C
    1 a1 b1 c1
    2 a2 b2 c2

    AB <- data.frame(A = c("a1", "a2"), B = c("b1", "b2"))
    BA <- data.frame(B = c("b1", "b2"), A = c("a1", "a2"))
    AC <- data.frame(A = c("a1", "a2"), C = c("c1", NA))
    CA <- data.frame(C = c("c1", NA), A = c("a1", "a2"))

    BC <- data.frame(B = c("b1", "b2"), C = c("c1", "c2"))
    CB <- data.frame(C = c("c1", "c2"), B = c("b1", "b2"))

    AB_ <- full_join(AB, BA)

    Joining, by = c("A", "B")

    AC_ <- full_join(AC, CA)

    Joining, by = c("A", "C")

    A <- full_join(AB_, AC_)

    Joining, by = "A"

    BC_ <- full_join(BC, CB)

    Joining, by = c("B", "C")

    BA_ <- full_join(BA, AB)

    Joining, by = c("B", "A")

    B <- full_join(BA_, BC_)

    Joining, by = "B"

    full_join(A, B)

    Joining, by = c("A", "B", "C")

       A  B    C
    1 a1 b1   c1
    2 a2 b2 <NA>
    3 a2 b2   c2
