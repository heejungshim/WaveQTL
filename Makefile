all:      
	cd src && $(MAKE)
	cd src && mv ./WaveQTL ../WaveQTL
	mkdir -p ./waveqtl2d/bin
	cp ./WaveQTL ./waveqtl2d/bin/WaveQTL
clean:
	cd src && $(MAKE) clean
	rm -rf ./waveqtl2d/bin
	rm -f ./WaveQTL
	
