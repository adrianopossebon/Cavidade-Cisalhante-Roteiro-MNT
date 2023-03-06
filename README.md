# Cavidade Cisalhante: Introdução

Aqui eu resolvo, numericamente, o problema do escoamento bidimensional de um fluido newtoniano e incompressível em uma cavidade cisalhante (tampa de cima se movendo). Trata-se de um problema clássico na Mecânica dos Fluidos Computacional. As equações da continuidade e de Navier-Stokes são resolvidas com um método de projeção explícito de primeira ordem no tempo. É utilizado o método das diferenças finitas juntamente com uma malha defasada. 

Os detalhes do problema e da implementação estão todos no roteiro

	MNT_Cavidade_Cisalhante.pdf


# Arquivos

O código está distribuído em 5 arquivos:

1. **cavity_main.py**

	Arquivo com a parte principal do código, incluindo a entrada dos parâmetros e o _loop_ principal.
	
1. **cavity_velocity.py**
	
	Arquivo com as funções relacionadas à velocidade. 
	
1. **cavity_pressure.py**

	Arquivo com as funções relacionadas ao cálculo da pressão. 

1. **cavity_general_functions.py**

	Arquivo com as funções auxiliares.

1. **cavity_plot.py**

	Arquivo com as funções usadas para plotar os resultados finais. 


# Setup

Versão do python e das bibliotecas que eu estou usando:

* python==3.8.5
* numpy==1.19.2
* pandas==1.1.3
* numba==0.51.2
* matplotlib==3.3.2


# Execução

Para rodar o código principal no Linux:
	
	python cavity_main.py

Também é possível rodar apenas o código que gera os gráficos, quando os arquivos de output já estão disponíveis. Neste caso, o comando é:

	python cavity_plot.py

Atenção: é preciso definir alguns parâmetros no arquivo cavity_plot.py antes de rodar. 