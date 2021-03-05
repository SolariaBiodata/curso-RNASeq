
---
permalink: /sesion04.html
---
![alt text](https://solariabiodata.com.mx/images/solaria_banner.png "Soluciones de Siguiente Generación")
# Curso de Análisis de Datos de RNA-Seq

## Sección 04: Mapeo de Referencia y Cuantificación

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
9. Finalmente, realizaremos los perfiles de cuantificación de cada muestra, con la noción que usaremos Ballgown para el proceso de expresión Diferencial.
    ~~~
    stringtie -e -B -p 24 -G stringtie_merged.gtf \
      -o ballgown/SampleA/transcripts.gtf sampleA.sort.bam
    ~~~
    *** NOTA: *** Si es deseable, podemos trabajar con una operación de ciclo FOR para todas las muestras de una carpeta.
    ~~~
    for i in $(ls *sort.bam); do
      stringtie -e -B -p 24 -G stringtie_merged.gtf \
      -o ballgown/${i/.sort.bam}/transcripts.gtf $i;
    done
    ~~~

## Ejercicio 02: Ensayo de Expresión Diferencial
### Descripción
Una vez que se han generado los perfiles de expresión, utilizaremos el lenguaje de programación R para generar la expresión diferencial.

### Instrucciones
1. Crearemos un archivo de metadatos denominado `pheno.csv`, con el contenido siguiente como ejemplo:
    ~~~
    ids,sample,experiment
    SampleA_S7_L001,SampleA,control
    SampleA_S7_L002,SampleA,control
    SampleA_S7_L003,SampleA,control
    SampleA_S7_L004,SampleA,control
    SampleB_S8_L001,SampleB,experimental
    SampleB_S8_L002,SampleB,experimental
    SampleB_S8_L003,SampleB,experimental
    SampleB_S8_L004,SampleB,experimental
    ~~~
2. Cargaremos la consola de R, o ejecutaremos RStudio en el caso de contar con interfaz gráfica.
3. En la consola de R o en la construcción del script de análisis, cargaremos las librerías necesarias.
    ~~~
    library(tidyverse)
    library(dplyr)
    library(genefilter)
    library(RSkittleBrewer)
    library(ballgown)
    ~~~
4. Cargamos los objetos del archivo pheno.csv y toda la carpeta de resultados que generamos en el punto 9 del ejercicio anterior
    ~~~
    pheno = read.csv("pheno.csv")
    bg = ballgown(dataDir = "ballgown",samplePattern ="Sample",pData = pheno)
    ~~~
5. Utilizaremos gexpr para los perfiles de expresión diferencial de este objeto:
    ~~~
    gexpr(bg)
    ~~~
6. Filtraremos los datos para FPKM mayores a 1
    ~~~
    bgf = subset(bg,"rowVars(gexpr(bg)) > 1", genomesubset=T)
    ~~~
7. Definiremos un intervalo de confidencia para la prueba estadística de expresión Diferencial
    ~~~
    alpha = 0.01
    ~~~
8. Ensayo de Expresion Diferencial
    ~~~
    res = stattest(bgf, feature = "gene", covariate = "experiment", getFC = T, meas = "FPKM")
    ~~~
9. Resultados de la expresión diferencial.
    ~~~
    res %>% head()
    ~~~
10. Agregamos las siguientes columnas.
    ~~~
    res$stat = ifelse(res$qval<alpha,yes=TRUE,no=FALSE) # Definimos alpha como valor de corte para asumir diferencia significativa
    res$log2fc = log(res$fc,base=2)                     # Agregamos la razón de cambio en log2
    res$log10qvalue = -log(res$qval,base=10)            # Agregamos el valor estadístico en -log10
    ~~~
11. Por último, filtraremos los datos con expresión significativa
    ~~~
    diff_genes = res %>% filter(stat == TRUE)
    diff_ids = diff_genes$id
    ~~~
