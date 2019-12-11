TAG = hurwitzlab/trim_galore:0.6.5

img:
	docker build --tag=$(TAG) .

run:
	docker run --rm $(TAG)

shell:
	docker run --rm -it $(TAG) bash

test:
	docker run --rm -v /data:/data $(TAG) run_centrifuge.py -q /data/kyclark/dolphin/test/Dolphin_1.fasta -o /data/kyclark/centrifuge-out -I /data/centrifuge-indexes

dolphin:
	docker run --rm -v /Users/kyclark/work/data:/data $(TAG) run_centrifuge.py -q /data/dolphin/fasta -o /data/dolphin/centrifuge-out -I /data/centrifuge-indexes -i p_compressed --procs 1

push: 
	docker push $(TAG)
