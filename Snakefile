from os import listdir
from os.path import isfile, join

configfile: "config.yaml"

file = [f for f in listdir(config["samples_dir"]) if isfile(join(config["samples_dir"], f))]
files = list(set([f.split('_')[0] for f in file]))
dir = config["samples_dir"]


rule all:
    input:
        "multiqc_reports/multiqc_report.html",
        "trimmed_multiqc_reports/multiqc_report.html",
        expand("krona/kraken/krona_{s}_1.html", s =files),
        expand("krona/kaiju/{s}_1.out.html", s=files),
        directory("mapDamage/results_mtdna/"),
        directory("prokka/mtdna/"),
        "haplogrep/haplogroups.txt",
        "samtools_stats/statistics_mtdna.txt",
        "phylo/results/phylo_tree.png"
        

"""
Diese Regel führt eine Qualitätskontrolle mittels FastQC auf fastq Dateien durch. Als Ausgabe wird für jedes fastq-File
eine zip und html mit den Zusammenfassungen im Ordner fastqc_reports abgelegt.
"""
rule fastqc:
    input:
        expand(config["samples_dir"] + "{sample_name}_{n}.fastq",sample_name=files, n=["1","2"])
    output:
        expand("fastqc_reports/{sample_name}_{n}_fastqc.{t}", sample_name=files, n=[1,2], t=["zip","html"])
    threads:
        config["threads"]
    shell:
        "fastqc {input} -t {threads} -o fastqc_reports/"


"""
Diese Regel fasst fastqc-Reports im fastqc_reports Ordner zusammen und legt das Ergebnis im Ordner multiqc_reports ab. 
"""
rule multiqc:
    input:
        expand("fastqc_reports/{sample_name}_{n}_fastqc.{t}", sample_name=files, n=[1,2], t=["zip","html"])
    output:
        "multiqc_reports/multiqc_report.html"
    params:
        multiqc_output="multiqc_reports/"
    shell:
        "multiqc -o {params.multiqc_output} {input}"

"""
In dieser Regel werden fastq-Files mit dem Tool trim_galore getrimmt. Dieses Trimming beinhaltet Qualitäts- und Adaptertrimming.
Die getrimmten Sequenzen und der Bericht werden im Ordner trimmed_seqs gespeichert.
"""
rule trimming:
    input:
        samples = expand(config["samples_dir"] + "{sample_name}_1.fastq", sample_name=files)
    output:
        expand("trimmed_seqs/{sample_name}_1_trimmed.fq", sample_name=files),
        expand("trimmed_seqs/{sample_name}_1.fastq_trimming_report.txt", sample_name=files),
    threads:
        config["threads"]
    shell:
        "trim_galore -j {threads} {input.samples} -o trimmed_seqs/"

"""
Mit dieser Regel wird FastQC auf die getrimmten Sequenzen ausgeführt.
Die Ausgabe wird in den fastqc_reports Ordner abgelegt.
"""
rule fastqc_trimmed:
    input:
        expand("trimmed_seqs/{sample_name}_1_trimmed.fq",sample_name=files)
    output:
        expand("fastqc_reports/{sample_name}_1_trimmed_fastqc.{t}", sample_name=files, t=["zip","html"])
    threads: 
        config["threads"]
    shell:
        "fastqc {input} -t {threads} -o fastqc_reports/"

"""
Diese Regel führt MultiQC auf den FastQC Reports der getrimmten Sequenzen aus.
Die Zusammenfassung wird im Ordner trimmed_multiqc_reports abgelegt.
"""
rule multiqc_trimmed:
    input:
        expand("fastqc_reports/{sample_name}_1_trimmed_fastqc.{t}", sample_name=files, t=["zip","html"])
    output:
        "trimmed_multiqc_reports/multiqc_report.html"
    params:
        multiqc_output="trimmed_multiqc_reports/"
    shell:
        "multiqc -o  {params.multiqc_output} {input}"


"""
Diese Regel führt das Tool Kraken2 aus. Hierbei werden die getrimmten Sequenzen und eine Kraken Datenbank als Eingabe verwendet.
Für jedes Sample wird eine Klassifikation, Zusammenfassung und Report im Kraken Ordner gespeichert.
"""
rule kraken2:
    input:
        data="trimmed_seqs/{sample}_trimmed.fq",
        db= config["kraken_db"]
    output:
        classification="kraken/classification_{sample}.txt",
        summary="kraken/summary_{sample}.txt",
        report="kraken/report_{sample}.txt"
    threads:
        config["threads"]
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.classification} --use-mpa-style {input.data} > {output.summary}) 2> {output.report}"


"""
Diese Regel konvertiert eine Zusammenfassung von Kraken2 in eine Krona Ausgabe im HTML Format.
Die Eingabe für diesen Aufruf ist die Zusammenfassung der einzelnen Samples.
"""
rule kraken2krona:
    input:
        "kraken/summary_{sample}.txt"
    output:
        "krona/kraken/krona_{sample}.html"
    shell:
        "ktImportTaxonomy -t 5 -m 3 -o {output} {input}"


"""
Diese Regel führt das Tool Kaiju aus. Als Eingabe werden die getrimmten Sequenzen und die in der Config festgelegte Kaiju Datenbank mitgegeben.
Die erste Ausgabe ist ein eine .out Datei von Kaiju. Anschließend wird das Tool kaiju2krona ausgeführt, welches die Ausgabe von Kaiju in ein Krona Format umwandelt.
Zuletzt wird ktImportText ausgeführt, welches den Krona Output in HTML Format umwandelt.
"""
rule kaiju:
    input:
        trimmed_fastq="trimmed_seqs/{sample}_trimmed.fq",
        database_nodes=config["kaiju_db"] + "/nodes.dmp",
        database_names=config["kaiju_db"] + "/names.dmp"
    output:
        kaiju_output="kaiju/{sample}.out",
        krona_output="krona/kaiju/{sample}.out.krona",
        html_output="krona/kaiju/{sample}.out.html"
    shell:
        "kaiju -z 12 -t {input.database_nodes} -f /vol/biodb/local_databases/MGX/kaiju/current/kaiju_db.fmi -i {input.trimmed_fastq} -v -o {output.kaiju_output};"
        "kaiju2krona -t {input.database_nodes} -n {input.database_names} -i {output.kaiju_output} -o {output.krona_output};"
        "ktImportText -o {output.html_output} {output.krona_output}"


"""
Diese Regel baut einen Index aus dem humanen Mitochondrial Referenzgenom. Hierbei wird eine .fa Datei des humanen Referenzgenoms als Eingabe verwendet.
Im Ordner bowtie2 wird ein weiterer Ordner namens index_mtdna angelegt, in welchem die Dateien abgelegt werden. 
"""
rule bowtie2_build_mtdna:
    input:
        "ref_genome/ref_genome_human_mitochondrial.fa"
    params:
        basename= "bowtie2/index_mtdna/index_mtdna"
    output:
        expand("bowtie2/index_mtdna/index_mtdna.{entry}.bt2", entry=["1", "2", "3", "4", "rev.1", "rev.2"])
    shell:
        "bowtie2-build {input} {params.basename}"              

"""
Diese Regel mappt die getrimmten Sequenzen auf das mitochondriale Referenzgenom. Dabei werden die gemappten und ungemappten Sequenzen im Dateiformat fastq gespeichert. Zudem wird eine Zusammenfassung des Mapping erstellt. Die gemappten Sequenzen werden außerdem als .sam Datei in bowtie2/mtdna_mapped abgelegt. Um Ressourcen zu sparen wird die .sam Datei nur temporär gespeichert und nach Beedingung des Workflows entfernt.
"""
rule bowtie2_map_mtdna:
    input:
        indices = expand("bowtie2/index_mtdna/index_mtdna.{entry}.bt2", entry=["1", "2", "3", "4", "rev.1", "rev.2"]),
        sample= "trimmed_seqs/{trimmed_sample}_1_trimmed.fq"
    output:
        mapped = "bowtie2/mtdna_mapped/{trimmed_sample}_1.fastq", 
        alignment = temp("bowtie2/mtdna_mapped/aligned_{trimmed_sample}_1.sam"),
        summary = "bowtie2/mtdna_summary/{trimmed_sample}_1.txt",
        unmapped = "bowtie2/mtdna_unmapped/{trimmed_sample}_1.fastq"
    params:
        index = "bowtie2/index_mtdna/index_mtdna"
    threads:
        config["threads"]
    shell:
        "(bowtie2 -p {threads} -x {params.index} -U {input.sample} --un {output.unmapped} --al {output.mapped} -S {output.alignment}) 2> {output.summary}"

"""
Diese Regel verwendet samtools und wandelt die .sam Dateien, aus bowtie2 in .bam Dateien um.
Diese Dateien werden im gleichen Ordner wie die .sam Dateien abgelegt.
"""
rule sam2bam_mtdna_map:
    input:
        "bowtie2/mtdna_mapped/aligned_{sample}.sam"
    output:
        "bowtie2/mtdna_mapped/aligned_{sample}.bam"
    shell:
        "samtools view -bh {input} > {output}"


"""
Diese Regel führt alle gemappten Sequenzen im .bam Format in einer Datei zusammen.
Diese Datei wird im gleichen Ordner abgelegt.
"""
rule merge_mtdna_bam:
    input:
        expand("bowtie2/mtdna_mapped/aligned_{s}_1.bam", s=files)
    output:
        "bowtie2/mtdna_mapped/merged_mtdna_mapped.bam"
    shell:
        "samtools merge -o {output} {input}"

"""
Diese Regel sortiert die zusammengeführte .bam Datei der mtDNA. 
Die sortierte .bam Datei wird im gleichen Ordner abgelegt wie die unsortierte Datei.
"""
rule bam_sort_mtdna_map:
    input:
        "bowtie2/mtdna_mapped/merged_mtdna_mapped.bam"
    output:
        "bowtie2/mtdna_mapped/sorted_merged_mtdna.bam"
    shell:
        "samtools sort -O bam -o {output} {input}"


"""
Diese Regel legt einen bam Index für eine sortierte .bam Datei an.
Dieser Index wird im gleichen Ordner abgelegt wie die .bam Datei.
"""
rule bam_index_mtdna_map:
    input:
        "bowtie2/mtdna_mapped/sorted_merged_mtdna.bam"
    output:
        "bowtie2/mtdna_mapped/sorted_merged_mtdna.bam.bai"
    shell:
        "samtools index {input} -o {output}"


"""
Diese Regel legt die Consensus Sequenz aus der zusammengeführten und sortierten .bam Datei an.
Die resultiernde Sequenz wird als .fa Datei im ordner consensus/mtdna abgelegt.
"""
rule consensus_mtdna:
    input:
        bam="bowtie2/mtdna_mapped/sorted_merged_mtdna.bam",
        bai="bowtie2/mtdna_mapped/sorted_merged_mtdna.bam.bai"
    output:
        "consensus/mtdna/merged_mtdna_consensus.fa"
    shell:
        "samtools consensus -a --show-ins no {input.bam} -o {output}"


"""
Diesie Regel führt eine Genomannotation mittels Prokka aus. Die Eingabedatei ist hierbei die Consensus Sequenz, welche in der vorherigen Regel erstellt wurde. Als Ausgabe wird ein Ordner angelegt, in welchem diverse Dateiformate wie ".gff" zu der Genomannotation aufzufinden sind.
Die übergebenen Parameter geben an, dass es sich um mitochondirale DNA handelt.
"""
rule annotation_mtdna:
    input:
        "consensus/mtdna/merged_mtdna_consensus.fa"
    output:
        directory("prokka/mtdna/")
    params:
        mode="Mitochondria",
        prefix="mtdna"
    shell:
        "prokka {input} --kingdom {params.mode} --prefix {params.prefix} --outdir {output} --force"


"""
Diese Regel führt das Tool mapDamage aus, welches Schäden von ancient DNA quantifiziert. Die Eingabe dafür ist ein Referenzgenom und eine .bam Datei mit allen Reads.
Als Ausgabe wird ein Ordner angelegt mit diversen Ausgaben wie Plots, Statistiken und Berichten.
"""
rule mapDamage_mtdna:
    input:
        bam = "bowtie2/mtdna_mapped/sorted_merged_mtdna.bam",
        ref = "ref_genome/ref_genome_human_mitochondrial.fa"
    output:
        folder = directory("mapDamage/results_mtdna")
    shell:
        "mapDamage -i {input.bam} -r {input.ref} -d {output.folder}"


"""
Diese Regel ermittelt die Statistiken der zusammengefügten .bam Datei.
Die daraus folgende Ausgabe wird in eine .txt Datei im Ordner samtools_stats abgelegt.
"""
rule samtools_stats:
    input:
        "bowtie2/mtdna_mapped/sorted_merged_mtdna.bam"
    output:
        "samtools_stats/statistics_mtdna.txt"
    shell:
        "samtools stats {input} > {output}"


"""
Diese Regel führt das Tool Haplogrep aus. Die Eingabe hierfür ist die Consensus Sequenz der mtDNA.
Als Ausgabe wird eine .txt Datei in den Ordner haplogrep abgelegt, welche die Haplogruppen vorhersagt.
"""
rule haplogrep:
    input:
        "consensus/mtdna/merged_mtdna_consensus.fa"
    output:
        "haplogrep/haplogroups.txt"
    shell:
        "haplogrep classify --in {input} --format fasta --out {output}"

"""
Diese Regel führt einige Genome aus der amtDB mit der Consensus Sequenz zusammen. Dabei wird die Ausgabe so prozessiert, dass Header und Sequenzen in der Datei richtig gesetzt werden.
"""
rule merge_consensus:
    input:
        multiple_sequences="phylo/sequences/prepared_test_genomes.fasta",
        consensus="consensus/mtdna/merged_mtdna_consensus.fa"
    output:
        merged_sequences="phylo/sequences/merged_phylo_sequences.fasta",
        prepared_sequences="phylo/sequences/prepared_genomes.fasta"
    shell:
        "cat {input.multiple_sequences} {input.consensus} > {output.merged_sequences};"
        "sed '/^$/d' {output.merged_sequences}  > {output.prepared_sequences}"


"""
Diese Regel führt das Tool Clustal Omega aus. Die Eingabe für dieses Tool sind mehrere zusammengeführte Sequenzen. 
Mit diesem Tool wird ein multiples Sequenzalignment durchgeführt, welches anschließend als .fasta Datei gespeichert wird.
"""
rule clustalo:
    input:
        "phylo/sequences/prepared_genomes.fasta"
    output:
        "phylo/sequences/alignment.fasta"
    shell:
        "clustalo -i {input} -o {output} -v"

"""
Diese Regel führt ein Python Skript aus, welches einen phylogenetischen Baum aus dem multiplen Sequenzalignment anlegt. Dafür wird die Ausgabe von Clustal Omega benötigt.
Das Skript erstellt einen phylogenetischen Baum aus dem Alignment im .png Format und speichert dieses im Ordner phylo/results/.
"""
rule phylo_tree:
    input:
        "phylo/sequences/alignment.fasta"
    output:
        "phylo/results/phylo_tree.png"
    shell:
        "python3 phylo/scripts/tree.py"



