---
permalink: /sesion02.html
---
![alt text](https://solariabiodata.com.mx/images/solaria_banner.png "Soluciones de Siguiente Generación")
# Curso de Análisis de Datos de RNA-Seq

## Sección 02: Revisión y Filtrado de Calidad

### Descripción
En esta sesión aprenderemos lo necesario para descargar secuencias de NCBI, filtrar las búsquedas con opciones avanzadas y empezar a realizar comparaciones en resultados de alineamientos.

### Requisitos

Para poder realizar este ejercicio, necesitaremos:

1. Datos de secuencias:
    - Puedes usar tus propias secuencias, en formato FASTQ o FASTQ.GZ
2. Acceso a los siguientes recursos de internet:
    - Página de [NCBI](https://www.ncbi.nlm.nih.gov/)
    - Página de [SRA](https://www.ncbi.nlm.nih.gov/sra)
3. Sofware Recomendable para esta sesión:
    - Terminal (Mac o Linux)
    - Putty (Windows)
    - WinSCP (Windows)

## Ejercicio 01: Remoción de Adaptadores
### Descripción
Dada la heterogeneidad de causas para una baja calidad, se deben tomar a consideración diferentes escenarios para la limpieza y filtrado de secuencias.

### Instrucciones (Parte 01)
1. Para remover adaptadores usando Trimmomatic:
    ~~~
    Trimmomatic SE biblioteca_1_R1.fastq.gz Trimm-biblioteca_1_R1.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10
    ~~~
### Instrucciones (Parte 02)
2. Removiendo adaptadores con la utilidad de fastx-toolkit:
    ~~~
    fastx_clipper -a ATTTGGTACGTA -i mislecturas.fastq -o mislecturas_noadapter.fastq
    ~~~

### Instrucciones (Parte 03)
3. Usando cutadapt para remover adaptadores
    ~~~
    cutadapt -a ATTTGGTACGTA -o mislecturas_noadapter.fastq mislecturas.fastq
    ~~~

## Ejercicio 02: Corte por Calidad
### Descripción
Dada la heterogeneidad de causas para una baja calidad, se deben tomar a consideración diferentes escenarios para la limpieza y filtrado de secuencias.

### Instrucciones (Parte 01)
1. Para remover por Calidad en el extremo 5' de la librería:
    ~~~
    Trimmomatic SE biblioteca_1_R1.fastq.gz Trimm-biblioteca_1_R1.fastq.gz LEADING:20
    ~~~

### Instrucciones (Parte 02)
2. Para remover por Calidad en el extremo 3' de la librería
    ~~~
    Trimmomatic SE biblioteca_1_R1.fastq.gz Trimm-biblioteca_1_R1.fastq.gz TRAILING:20
    ~~~

### Instrucciones (Parte 03)
3. Removiendo adaptadores con la utilidad de fastx-toolkit:
    ~~~
    fastq_quality_filter -q 30 -p 90 -i mislecturas.fastq -o mislecturas_noadapter.fastq
    ~~~

## Ejercicio 02: Corte Fijo
### Descripción
Si por lo general, los extremos de los datos de NGS tienden a bajar la calidad, podrían recortarse los tamaños de lectura.

### Instrucciones (Parte 01)
1. Para cortar 20 bases en el extremo 5' de la librería:
    ~~~
    Trimmomatic SE mislecturas_R1.fastq.gz headcrop_mislecturas_R1.fastq.gz HEADCROP:20
    ~~~

### Instrucciones (Parte 02)
2. Para cortar 20 bases en el extremo 3' de la librería
    ~~~
    Trimmomatic SE mislecturas_R1.fastq.gz crop_mislecturas_R1.fastq.gz CROP:20
    ~~~

### Instrucciones (Parte 03)
3. Cortando 5 bases con la utilidad de fastx-toolkit:
    ~~~
    cutadapt -u 5 -o mislecturas.fastq cut5_mislecturas.fastq
    ~~~

## Ejercicio 03: Ventana Deslizante
### Descripción
Esta estrategia nos permite minimizar la pérdida de información de buena calidad

### Instrucciones
1. Para realizar ventanas deslizantes de 20 bases con un mínimo de Q30:
    ~~~
    Trimmomatic SE mislecturas.fastq.gz sw_mislecturas_R1.fastq.gz SLIDINGWINDOW:20:30
    ~~~
2. Para realizar ventanas deslizantes de 20 bases con un mínimo de Q30 y deslizamiento cada 5 bases:
    ~~~
    Trimmomatic SE mislecturas.fastq.gz sw_mislecturas_R1.fastq.gz SLIDINGWINDOW:20:30:5
    ~~~

## Ejercicio 04: Modos Adicionales de Trimmomatic
### Descripción
Algunas funcionalidades adicionales que pueden ser muy útiles en Trimmomatic son el uso de librerías pareadas y Concatenación de distintas opciones de filtro:
### Instrucciones
1. Para remover lecturas menores a 150 bases:
    ~~~
    Trimmomatic SE mislecturas_R1.fastq.gz minlen_mislecturas_R1.fastq.gz MINLEN:150
    ~~~
2. Para usar múltiples estrategias de filtrado en una sola instrucción:
    ~~~
    Trimmomatic SE mislecturas_R1.fastq.gz filter_mislecturas_R1.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 \
    SLIDINGWINDOW:4:25 MINLEN:150

    ~~~
3. Para filtrado de lecturas pareadas:
    ~~~
    Trimmomatic PE mislecturas_R1.fastq.gz mislecturas_R2.fastq.gz \
    R1_paired.fq.gz R1_unpaired.fq.gz R2_paired.fq.gz R2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 \
    SLIDINGWINDOW:4:25 MINLEN:150
    ~~~
 

## Ejercicio 05: Ejecución serial de Trimmomatic con todos los filtros para todas las muestras:
### Descripción
Usando las lecturas propias, podemos seleccionas valores específicos para `ILLUMINACLIP`, `LEADING`, `TRAILING`, `SLIDINGWINDOW` y `MINLEN` de tal forma que generemos 

### Instrucciones

1. En general tenemos que realizar este flujo para cada pareja de lecturas en todas las muestras:
    ~~~
    cd raw
    
    Trimmomatic PE muestra_R1.fastq.gz muestra_R2.fastq.gz \
    muestra_1P.fq muestra_1U.fq muestra_2P.fq muestra_2U.fq \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:5:25 MINLEN:60 
    ~~~

2. Una vez que hemos generado todas las lecturas limpias, simplemente movemos los datos a una carpeta diferente con las lecturas limpias:
    ~~~
    cd ..
    mkdir libs
    mv raw/*P.fq raw/*U.fq libs/
    ~~~


