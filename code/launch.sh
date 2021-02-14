#!/bin/bash
#launch.sh

echo "Ejecutando flujo de trabajo"
if [ -d "software/deps_r" ];then
	echo "Introduce la ruta donde se localiza el proyecto"
	read filename
	# echo "Ejecutando archivo de configuración"
	# Rscript configure.R 
	echo "Ejecutando archivo principal"
	Rscript human_covid_network.R $filename
else
    echo "Instalando librerías"
    sudo bash setup.sh
    echo "Introduce la ruta donde se localiza el proyecto"
    read filename
    echo "Ejecutando archivo principal"
    Rscript human_covid_network.R $filename
fi
echo "Ejecucion finalizada"