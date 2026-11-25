version 1.3

## Count how many DNase I hypersensitivity peaks overlap features of interest.
workflow peak_overlaps {
    meta {
        description: "Counts peaks overlapping a set of genomic features, one sample at a time"
        outputs: {
            summary: "A tab-separated table with one row per sample",
        }
    }

    parameter_meta {
        features: "BED file of features to intersect with (e.g. exons)"
        samples: "BED files of DNase I hypersensitivity peaks, one per sample"
    }

    input {
        File features
        Array[File] samples
    }

    scatter (sample in samples) {
        call count_overlaps {
            features,
            peaks = sample,
        }
    }

    call summarize {
        counts = count_overlaps.counts,
    }

    output {
        File summary = summarize.table
    }
}

task count_overlaps {
    meta {
        description: "Counts the peaks in one sample that overlap at least one feature"
        outputs: {
            overlapping: "Peaks that overlap at least one feature",
            counts: "One line: sample name, total peaks, overlapping peaks",
        }
    }

    parameter_meta {
        features: "BED file of features"
        peaks: "BED file of peaks for one sample"
    }

    input {
        File features
        File peaks
    }

    String name = basename(peaks, ".bed")

    command <<<
        set -euo pipefail
        bedtools intersect -u -a "~{peaks}" -b "~{features}" > overlapping.bed
        printf '%s\t%s\t%s\n' \
            "~{name}" \
            "$(wc -l < "~{peaks}")" \
            "$(wc -l < overlapping.bed)" \
            > counts.tsv
    >>>

    output {
        File overlapping = "overlapping.bed"
        File counts = "counts.tsv"
    }

    requirements {
        container: "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
        cpu: 1
        memory: "64 MiB"
    }
}

task summarize {
    meta {
        description: "Combines per-sample counts into a single table"
        outputs: {
            table: "Tab-separated summary with a header row",
        }
    }

    parameter_meta {
        counts: "Per-sample count files produced by `count_overlaps`"
    }

    input {
        Array[File] counts
    }

    command <<<
        set -euo pipefail
        printf 'sample\ttotal_peaks\toverlapping_peaks\n' > summary.tsv
        cat ~{sep(" ", counts)} >> summary.tsv
    >>>

    output {
        File table = "summary.tsv"
    }

    requirements {
        container: "ubuntu:24.04"
        cpu: 1
        memory: "64 MiB"
    }
}
