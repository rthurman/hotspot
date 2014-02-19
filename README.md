hotspot
=======

Hotspot is a program for identifying regions of local enrichment of
short-read sequence tags mapped to the genome using a binomial
distribution model. Regions flagged by the algorithm are called
"hotspots." The algorithm utilizes a local background model that
automatically normalizes for large regions of elevated tag levels due
to, for example, copy number effects. Hotpsot is otherwise able to
detect regions of enrichment of highly-variable size, making it
applicable to both broad and highly-punctate signals. We have applied
it extensively to DNase-seq and ChIP-seq data, including transcription
factor (CTCF) and histone modification (H3K4me3, H3K36me3, H3K27me3)
data.  This distribution also includes scripts for computing SPOT
(Signal Portion of Tags), a quality measure for short-read sequence
experiments. SPOT is simply the percentage of all tags that fall in
hotspots.

Further documentation can be found in the distribution documentation,
and on the Hotspot/SPOT web site,

[http://www.uwencode.org/proj/hotspot](http://www.uwencode.org/proj/hotspot)


