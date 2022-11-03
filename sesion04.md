
---
permalink: /sesion04.html
---
![alt text](https://solariabiodata.com.mx/images/solaria_banner.png "Soluciones de Siguiente Generación")
# Curso de Análisis de Datos de RNA-Seq

## Sección 04: Mapeo de Referencia, Cuantificación y Expresión Diferencial

### Descripción
En esta sesión, realizaremos el perfil de expresión génica a partir de un transcriptoma de referencia (de la sección anterior).
### Requisitos

Para poder realizar este ejercicio, necesitaremos:

1. Datos de secuencias:
    - Puedes usar tus propias secuencias, en formato FASTQ o FASTQ.GZ datos limpios
2. Sofware Recomendable para esta sesión:
    - Terminal (Mac o Linux)
    - Putty (Windows)
    - WinSCP (Windows)

## Ejercicio 01: Mapeo con una Referencia
### Descripción
Es necesario concatenar todas las lecturas para hacer un ensamble completo. Una vez hecho eso se puede usar Trinity para hacer el ensamble.

### Instrucciones
1. Primero copiaremos el archivo de Contigs generados por Trinity, en una carpeta arriba de las librerías de secuenciación, sobre el espacio de trabajo.
    ~~~
    cp trinity_out/Trinity.fasta .
    ~~~
2. Prepararemos el índice de Burrows Wheeler para mapear transcriptomas.
    ~~~
    hisat2-build Trinity.fasta Trinity
    ~~~
3. Crearemos una carpeta donde depositaremos los mapeos de referencia
    ~~~
    mkdir mapping
    ~~~
4. Después de esto, prepararemos el mapeo con el comando hisat2.
    ~~~
    hisat2 --dta -x Trinity -1 libs/sampleA_1P.fq -2 libs/sampleA_2P.fq -S mapping/sampleA.sam
    ~~~
    *** NOTA: *** Si es deseable, podemos trabajar con una operación en lote para todas las muestras de una carpeta.
    ~~~
    for i in  $(ls libs/*1P.fq); do \
      hisat2 --dta -x Trinity \
      -1 $i -2 ${i/1P/2P} \
      -S ${i/_1P.fq}
    ~~~
    En este caso, si se ejecuta tal cual, hay que mover los mapeos que residen en `libs/` y hay que copiarlos en `mapping/`
    ~~~
    mv libs/*sam mapping/
    ~~~
    Ahora, cambiaremos el espacio de trabajo hacia la carpeta de mapeos.
    ~~~
    cd mapping
    ~~~
5. Para realizar el Ordenado e Indexado, primero hay que convertir el archivo de SAM a BAM.
    ~~~
    samtools view -bS -T Trinity.fasta sampleA.sam > sampleA.bam
    ~~~
    Realizaremos el ordenamiento
    ~~~
    samtools sort sampleA.bam > sampleA.sort.bam
    ~~~
    Y finalmente el Indexado
    ~~~
    samtools index sampleA.sort.bam
    ~~~
    *** NOTA: *** Si es deseable, podemos trabajar con una operación de ciclo FOR para todas las muestras de una carpeta.
6. Cuantificación de la Abundancia: Realizaremos la cuantificación de los transcritos mapeados en la referencia
    ~~~
    stringtie sampleA.sort.bam -p 16 -o sampleA.gtf
    ~~~
    *** NOTA: *** Si es deseable, podemos trabajar con una operación de ciclo FOR para todas las muestras de una carpeta.
    ~~~
    for i in $(ls *sort.bam); do
        stringtie $i -p 16 -o $i_out; done
    ~~~
7. Conjunción de todos los transcritos que se hayan cuantificado.
    ~~~
    ls *.gtf > gtf_todos.txt
    ~~~
8. Merging de Cuantificaciones, a partir de la lista guardada de los archivos GFF.
    ~~~
    stringtie --merge -o stringtie_merged.gtf gtf_todos.txt
    ~~~
9. Necesitamos crear una carpeta para depositar los resultados
    ~~~
    mkdir expression
    ~~~
10. Realizaremos los perfiles de cuantificación de cada muestra, con la noción que usaremos DESeq2 para el proceso de expresión Diferencial.
    ~~~
    stringtie -e -G stringtie_merged.gtf \
      -o expression/SampleA/transcripts.gtf sampleA.sort.bam
    ~~~
    *** NOTA: *** Si es deseable, podemos trabajar con una operación de ciclo FOR para todas las muestras de una carpeta.
    ~~~
    for i in $(ls *sort.bam); do
      stringtie -e -G stringtie_merged.gtf \
      -o expression/${i/.sort.bam}/transcripts.gtf $i;
    done
    ~~~
11. Por último, realizaremos el paso de generación de dos archivos: gene_count_matrix.csv y transcript_count_matrix.csv con PrepDE.py
    ~~~
    prepDE.py -i expression/
    ~~~

## Ejercicio 02: Ensayo de Expresión Diferencial
### Descripción
Una vez que se han generado los perfiles de expresión, utilizaremos el lenguaje de programación R para generar la expresión diferencial. Recuerda que puedes usar https://colab.to/R para crear una libreta de R, a la que debemos instalar el paquete DESeq2

### Instrucciones
0. Hay que asegurarnos que al crear una libreta de R; esto se hace con
    ~~~
    sessionInfo()
    ~~~
1. Primero, instalamos DESeq2 en la libreta de trabajo
    ~~~
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("DESeq2")
    ~~~
2. Crearemos un archivo de metadatos denominado `pheno.csv`, con el contenido siguiente como ejemplo:
    ~~~
    ids,sample,condition
    SampleA_S7_L001,SampleA,control
    SampleA_S7_L002,SampleA,control
    SampleA_S7_L003,SampleA,control
    SampleA_S7_L004,SampleA,control
    SampleB_S8_L001,SampleB,experimental
    SampleB_S8_L002,SampleB,experimental
    SampleB_S8_L003,SampleB,experimental
    SampleB_S8_L004,SampleB,experimental
    ~~~
3. En la celda de colab o en la construcción del script de análisis, cargaremos las librerías necesarias.
    ~~~
    library(tidyverse)
    library(dplyr)
    library(genefilter)
    library(RSkittleBrewer)
    library(DESeq2)
    ~~~
4. Cargamos los objetos del archivo pheno.csv y el archivo de gene_count_matrix.csv (o transcript_count_matrix.csv) del punto 11 del ejercicio anterior
    ~~~
    pheno = read.csv("pheno.csv")
    datos = read.table("/content/tabla_conteos.txt", header=T, sep="\t", row.names=1)
    ~~~
5. Utilizaremos DESeqDataSetFromMatrix() para los perfiles de expresión diferencial de este objeto:
    ~~~
    dds = DESeqDataSetFromMatrix(countData = as.matrix(datos), colData = pheno, design = ~ condition)
    ~~~
6. Para generar la expresión diferencial, invocamos la función principal DESeq()
    ~~~
    diffexp = DESeq(dds)
    res = results(diffexp)
    ~~~
7. Definiremos un intervalo de confidencia para la prueba estadística de expresión Diferencial
    ~~~
    alpha = 0.01
    ~~~
8. Ahora podemos realizar algunos filtros para extraer solo los genes expresados con significancia estadística
    ~~~
    stat_res = (res$pvalue<=0.01 & abs(res$log2FoldChange) >=2)
    ~~~
9. Resultados de la expresión diferencial.
    ~~~
    res %>% head()
    ~~~
10. Agregamos las siguientes columnas.
    ~~~
    res$stat = ifelse(res$qval<alpha,yes=TRUE,no=FALSE) # Definimos alpha como valor de corte para asumir diferencia significativa
    res$fc = 2^res$log2FoldChange                       # Agregamos la razón de cambio 
    res$log10fc = log(res$fc,base=10)                   # Agregamos la razón de cambio en log10
    res$log10qvalue = -log(res$qval,base=10)            # Agregamos el valor estadístico en -log10
    ~~~
11. Por último, filtraremos los datos con expresión significativa
    ~~~
    diff_genes = res %>% filter(stat == TRUE)
    diff_ids = diff_genes$id
    ~~~
12. Algunos gráficos útiles, pueden generarse como el Volcano Plot
    ~~~
    with(res, plot(log2FoldChange, -log10(pvalue), col="blue", pch=20, main="Volcano plot DESeq2"))
    ~~~
13. Y Otro gráfico es el MA Plot
    ~~~
    plotMA(res, main="MA-plot DESeq2")
    ~~~
