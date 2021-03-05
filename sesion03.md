---
permalink: /sesion04.html
---
![alt text](https://solariabiodata.com.mx/images/solaria_banner.png "Soluciones de Siguiente Generación")
# Curso de Análisis de Datos de RNA-Seq

## Sección 04: Perfiles de Expresión Diferencial

### Descripción
En esta sección aprenderemos a utilizar las lecturas limpias para realizar un ensamble de transcriptoma.

### Requisitos

Para poder realizar este ejercicio, necesitaremos:

1. Datos de secuencias:
    - Puedes usar tus propias secuencias, en formato FASTQ o FASTQ.GZ datos limpios
2. Sofware Recomendable para esta sesión:
    - Terminal (Mac o Linux)
    - Putty (Windows)
    - WinSCP (Windows)

## Ejercicio 01: Ensamble ***de novo***
### Descripción
Es necesario concatenar todas las lecturas para hacer un ensamble completo. Una vez hecho eso se puede usar Trinity para hacer el ensamble.

### Instrucciones (Parte 01)
1. Primero concatenamos las lecturas:
    ~~~
    cat libs/*1P.fq > all_1P.fq
    cat libs/*2P.fq > all_2P.fq
    ~~~
### Instrucciones (Parte 02)
2. Ejecutamos Trinity en modo paired end:
    ~~~
    Trinity --seqType fq --left all_1P.fq --right all_2P.fq --max_memory 84G --CPU 16 --output trinity_out
    ~~~

## Ejercicio 02: Revisión de la calidad del ensamble
### Descripción
Una vez que se han generado los contigs es importante revisar algunas métricas de calidad del ensamble.

### Instrucciones (Parte 01)
1. Revisiones generales, conteo de contigs, y revisando el tamaño de los contigs más largos:
    ~~~
    cd trinity_out/
    grep -c ">" Trinity.fasta
    grep ">" Trinity.fasta | cut -f2 -d" " | cut -f2 -d"=" | sort -h | tail
    cd ..
    ~~~

### Instrucciones (Parte 02)
2. Ejecución de quast:
    ~~~
    quast.py trinity_out/Trinity.fasta
    ~~~

En este punto, es importante notar que se ha generado una carpeta `quast_results`. Podemos copiar esa carpeta a nuestras computadoras locales para abrir los archivos html `results.html` con un explorador web.
