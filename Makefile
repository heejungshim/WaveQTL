all:      
	cd src && $(MAKE)
	cd src && mv ./WaveQTL ../WaveQTL
clean:
	cd src && $(MAKE) clean
	
