---
permalink: /sesion04.html
---
![alt text](https://solariabiodata.com.mx/images/solaria_banner.png "Soluciones de Siguiente Generación")
# Curso de Análisis de Datos de RNA-Seq

## Sesión Práctica del Día 06

### Descripción
En esta sesión sesión concluiremos con la miscelánea de visualización, utilizando librerías disponibles en los repositorios oficiales de R y BioconductoR

### Requisitos

Para poder realizar este ejercicio, necesitaremos:

1. Datos de secuencias:
    - Puedes usar tus propias secuencias, en formato FASTQ o FASTQ.GZ
2. Acceso a los siguientes recursos de internet:
    - Página de [NCBI](https://www.ncbi.nlm.nih.gov/)
    - Página de [SRA](https://www.ncbi.nlm.nih.gov/sra)
    - Página de [Ving](http://vm-gb.curie.fr/ving/)
3. Sofware Recomendable para esta sesión:
    - Terminal (Mac o Linux)
    - Putty (Windows)
    - WinSCP (Windows)

## Ejercicio 01: Visualización de Mapeos
### Descripción
A continuación, utilizaremos un script de construcción de imagenes de visualización, listas para publicación

### Instrucciones
1. Para ejecutar Ving.R sobre nuestros datos, usaremos la instrucción
    ~~~
    ving.R -o figure.png -c chromosome_name -S $start -E $end \
    -i -l -a anotacion.gff -r mRNA,CDS_parts -s line,box -C black,blue \
    1.mapped.bam n.mapped.bam
    ~~~

## Ejercicio 02: Generación de Gráficos en R
### Descripción
El siguiente ejercicio, usa la consola 'vanilla' de R, es decir, la implementación básica del lenguaje en la consola de comandos. Puedes seguir los siguientes ejercicios en una IDE, como lo es R Studio, en la sección de consola.

### Instrucciones (Parte 01)
1. Accederemos a la interfaz de R
    ~~~
    R
    ~~~
2. Una vez en la consola de R, descargaremos la siguiente fuente
    ~~~
    source("http://bioconductor.org/biocLite.R")
    ~~~
3. Instalación de la librería de cummeRbund
    ~~~
    biocLite(“cummeRbund”)
    ~~~
4. Cargado de la librería
    ~~~
    library(cummeRbund)
    ~~~
5. Posteriormente, crearemos un objeto y cargaremos la carpeta de resultados de cuffdiff
    ~~~
    cuff = readCufflinks(“diff_output”)
    ~~~
6. Visualizaremos la estructura del objeto cuff
    ~~~
    str(cuff)
    ~~~
### Instrucciones (Parte 02)
7. **Gráfico de Dispersión**: Representación de la posición de FPKM vs dispersión binomial por todas las réplicas biológicas
    ~~~
    disp = dispersionPlot(genes(cuff))
    ~~~
8. **Gráfico de Variabilidad**: Representación de la distribución de abundancia de transcritos vs CV^2
    ~~~
    scv = fpkmSCVPlot(genes(cuff))
    ~~~
9. **Gráfico de Distribución**: Representación de la distribución de FPKM por fracción de genes.
    ~~~
    dens = csDensity(genes(cuff))
    ~~~
10. **Gráfico de Caja y Bigotes**: Representación de la distribución de FPKM por total de condiciones.
    ~~~
    boxp = csBoxplot(genes(cuff))
    ~~~
11. **Gráfico de Matriz de Dispersión**: Representación de la distribución de FPKM entre las distintas condiciones evaluadas
    ~~~
    scatm = csScatterMatrix(genes(cuff))
    ~~~
12. **Gráfico de Dendograma de Outliers**: Representación de la relación y similitud entre el patrón de expresión de condiciones y réplicas.
    ~~~
    dende = csDendro(genes(cuff),replicates=T)
    ~~~
13. **Gráfico MA**: Representación del sesgo sistemático relacionado con la profundidad de secuenciación
    ~~~
    maplot = MAplot(genes(cuff),"Condicion1","Condicion2")
    ~~~
14. **Diagrama de Volcan**: Representación de la expresión diferencial de FPKM contra significancia estadística.
    ~~~
    volcamx = csVolcanoMatrix(genes(cuff))
    ~~~
15. **Mapa de Calor**: Representación de la expresión diferencial entre condiciones respecto a una escala de color
    ~~~
    data(sampleData)
    ~~~
    ~~~
    gene_names = getGenes(cuff,sampleIDs)
    ~~~
    ~~~
    heat<-csHeatmap(gene_names,cluster='both')
    ~~~
16. **Diagrama de Histograma de Expresión**: Representación de la expresión diferencial entre condiciones respecto a una escala de magnitud.
    ~~~
    barplot<-expressionBarplot(gene_names)
    ~~~
17. **Reducción Dimensional**: Representación de la asociación entre condiciones en un diagrama de componentes principales.
    ~~~
    rd_MDS = MDSplot(genes(cuff))
    ~~~
    ~~~
    rd_PCA = PCAplot(genes(cuff),"PC1","PC2")
    ~~~
### Instrucciones (Parte 03)
18. Exportación de Gráficos, primero escribiremos la función para exportar en formato PDF
    ~~~
    pdf(file=”mis_graficos.pdf”)
    ~~~
19. Escribiremos el nombre de todos los objetos que hemos creado
    ~~~
    disp; scv; dens; boxp; scatm; dende; maplot; volcam; heat; barplot; rd_MDS; rd_PCA
    ~~~
20. Terminaremos la exportación con lo siguiente
    ~~~
    dev.off()
    ~~~

## Ejercicio 03: Generación Alternativa de Heatmaps
### Descripción
Esta es una manera alternativa de crear un mapa de calor, usando las instrucciones nativas de R
### Instrucciones
1. Primero, en la terminal de Shell (Linux), usamos el siguiente script sobre el archivo genes_exp.diff
    ~~~
    cuffdiff_2_edgeR.pl --type "fpkm" --filename diff_output/genes.read_group_tracking > cuffdiff_genes.counts.matrix
    ~~~
2. Ahora, en la terminal de R, cargamos el archivo generado anteriormente
    ~~~
    cuffmx = read.table("cuffdiff_genes.counts.matrix",sep="\t",header=TRUE)
    ~~~
3. Guardaremos los nombres de los genes en un archivo aparte
    ~~~
    gene_names = as.vector(mx[,1])
    ~~~
4. Posteriormente, aplicaremos la instrucciónmx[,1]
    ~~~
    cuffmx = as.matrix(as.numeric(cuffmx[,-1]))
    ~~~
5. Por último, el heatmap se construye con
    ~~~
    heatmap_alt = heatmap(cuffmx,legend=gene_names)
    ~~~
6. Recuerda que se pueden exportar estos gráficos con las instrucciones del ejercicio 2, parte 3
