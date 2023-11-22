# Oetzi
Dieses Repository beinhaltet die Implementation eines Workflows, der eine Sequenzanalyse von Ötzi durchführt. Ötzi ist eine Gletschermumie, die im Jahr 1991 im österreichischen Ötztal entdeckt wurde. 
&nbsp;

* [Hintergrund](#hintergrund)

* [Aufgabenstellung](#aufgabenstellung)

* [Installation](#install)

* [Ausführen des Workflows](#workflow)

* [Ausgabe Workflow](#ergebnisse1)

* [Ergebnisse](#ergebnisse2)


### <a name="hintergrund"></a> Hintergrund
_Keller et al. (New insights into the Tyrolean Iceman's origin and phenotype as inferred by whole-genome sequencing)_ zufolge, welche das Genom von Ötzi im Jahr 2012 sequenziert haben, sollen unter anderem Sequenzen generiert worden sein, die Borrelia burgdoferi zugeordnet werden konnten. Durch diese Entdeckung entstand die Vermutung, dass Ötzi an Borreliose erkrankt sei. 

<img src="https://static.archaeologie-online.de/fileadmin/_processed_/6/9/csm_Oetzi_1abe4dab69.jpg"  width="320" height="300">

<img src="https://static.archaeologie-online.de/fileadmin/_processed_/e/2/csm_Oetzi-Portrait_2023_2f587cd7e2.jpg"  width="320" height="300">


Die oben dargestellten Abbildungen zeigen zum einem die gefundene Gletchermumie (links) [Bildquelle](https://static.archaeologie-online.de/fileadmin/_processed_/6/9/csm_Oetzi_1abe4dab69.jpg) sowie eine Illustration wie Ötzi, möglicherweise ausgesehen haben könnte (rechts) [Bildquelle](https://static.archaeologie-online.de/fileadmin/_processed_/e/2/csm_Oetzi-Portrait_2023_2f587cd7e2.jpg).



### <a name="aufgabenstellung"></a> Aufgabenstellung
Zur Untersuchung, ob Ötzi an Borreliose erkrankt war, soll im Rahmen dieses Projekts ein Workflow implementiert werden, der im Wesentlichen aus den folgenden vier Aufgabenstellungen besteht:
* Identifikation bakterieller Sequenzen in den Datensätzen
* Taxonomische Profilierung
* Rekonstruktion des Borrellia burgdorferi (Draft)-Genoms
* Rekonstruktion der mtDNA von Ötzi

Die Rekonstruktion des Borrellia burgdorferi (Draft)-Genoms sowie der mtDNA beinhaltet folgende Aufgaben:

* Mapping gegen das entsprechende Referenzgenom 
* Consensus Calling
* Annotation der erzeugten Contigs
* Bestimmung von 5’/3’ Substitutionsraten
* Phylogenetische Charakterisierung

Die hierfür verwendeten Daten sind im _European Nucleotide Archive_ (ENA) unter diesem [Link](https://www.ebi.ac.uk/ena/browser/view/PRJEB2830) zu finden.

Für die Rekonstruktion der Genome wurden zum einem das humane mitochondriale Referenzgenom aus der _National Center for Biotechnology Information_ ([NCBI](https://www.ncbi.nlm.nih.gov/nuccore/251831106)) sowie das von Keller et al. publizierte Burrelia burgdoferi Referenzgemom, das ebenfalls in der [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_011728.1/) hinterlegt ist, verwendet.


### <a name="install"></a> Installation
Das Repository mit dem von uns implementierten Workflows kann sowohl über `SSH` als auch über `HTTPS` heruntergeladen werden.
Eine Installation über `SSH` kann über folgendem Befehl ausgeführt werden:
```
git clone git@git.computational.bio.uni-giessen.de:rserradj/oetzi.git
```
Falls eine Installation über `HTTPS` erfolgen soll, kann dies mit dem nachfolgenden Befehl erreicht werden:
```
git clone https://git.computational.bio.uni-giessen.de/rserradj/oetzi.git
```
Die Installation des conda Enviroments setzt eine lokale Installation
von `conda` und `mamba` vorraus. Eine entsprechende [conda Installationsanleitung](https://conda.io/projects/conda/en/stable/user-guide/install/download.html) ist über den hier angegebenen Link zu finden. Die [mabma Installationsanleitung]() ist auf der folgenden Seite dokumentiert.

Nachdem `conda` und `mamba` erfolgreich installiert wurden, kann mittels `mamba` das von uns bereitgestellte `conda` Enviroment über die Kommandozeile installiert werden. Es empfiehlt sich, dies auf den `SLURM Cluster` auszuführen, da die Installation einige Zeit in Anspruch nehmen kann. Hierzu ist allerding eine Internetverbindung nötig.
``` 
mamba env create --file enviroment.yaml
```
Das conda Enviroment beinhaltet alle Pakete, die für den Workflow und Analyse benötigt werden und kann nach Installation über die Kommandozeile mit folgendem Befehl aktiviert werden:
```
conda activate oetzi
```
Nach der Aktivierung des conda Enviroments und vor Ausführung des Workflows muss folgender Befehl in der Kommandozeile ausgeführt werden.
```
ktUpdateTaxonomy.sh
```

### <a name="workflow"></a> Ausführen des Workflows
Der implementierte Workflow kann nach der Aktivierung des conda Enviroments über die Kommandozeile mit folgendem Befehl gestartet werden:
```
snakemake -j {Anzahl der Threads}
```
Das Ausführen des Workflows wird über den `SLURM Cluster` empfohlen.

Darüber hinaus kann der Nutzer über eine Konfigurationsdatei `config.yaml` den Pfad zu den Daten, Referenzgenomen und den Kraken2- sowie Kaiju-Datenbanken für die taxonimische Klassifikation angeben.
```
---
samples_dir: Pfad zu den Daten
burgdorferi_ref: Pfad zum Borrelia burgdoferi Referenzgenom
mtdna_ref: Pfad zum mtDNA Referenzgenom
kraken_db: Pfad zur kraken2 Datenbanken
kaiju_db: Pfad zur kajiu Datenbank
threads: Anzahl der Threads

```
### <a name="ergebnisse1"></a> Ausgabe Workflow
Nachdem der Workflow erfolgreich durchgelaufen ist, entsteht im Working Directory folgende Ordnerstruktur, die alle vom Workflow produzierten Ergebnisse beinhaltet:

```
Snakefile
├── fastqc_reeports
├── multiqc_reports
├── trimmed_seqs
├── trimmed_multiqc_reports
├── kraken
├── kaiju
├── krona
│   ├── kraken
│   └── kaiju  
├── bowtie2
│   ├── mtdna
│   │   ├── mapped
│   │   ├── unmapped
│   │   └── summary  
│   ├── index_burgdorferi 
│   └── index_mtdna 
├── consensus
│   └── mtdna
│       └── merged_mtdna_consensus.fa     
├── prokka
│   └── mtdna
├── mapDamage
│   └── reseults_mtdna
├── samtools_stats
│   └── statistics_mtdna.txt                     
└── phylo
    ├── sequences
    │   └── alignment.fasta
    └── reseults
        └── phylo_tree.png

```

### <a name="ergebnisse2"></a> Ergebnisse
Die visuelle Darstellung der Ergebnisse erfolgte mittels `HTML`. Hierfür muss die Datei `index.html` im lokalen Browser geöffnet werden.

Den dort dargestellten Ergebnisse ist zu entnehmen, dass Ötzi nicht an Borreliose erkrankt war, da aus den Daten kein Borrelia burgdoferi (Draft)-Genom rekonstruiert werden konnte.

Daüber hinaus konnte über die Bestimmung der Haplogruppe der mtDNA von Ötzi herausgefunden werden, dass dieser aus dem nördlichen Italien bzw. der Alpenregion stammte.

Eine genauere Darstellung aller Ergebnisse unseres Workflows ist in der Datei `HDDA.pdf` aufbereitet und dokumentiert.
