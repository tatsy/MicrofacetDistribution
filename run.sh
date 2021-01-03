#/bin/sh

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib beckmann \
    --alphax 0.1 \
    --alphay 0.1 \
    --output beckmann_walter07_0.1_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib beckmann \
    --alphax 0.5 \
    --alphay 0.1 \
    --output beckmann_walter07_0.5_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib beckmann \
    --alphax 0.5 \
    --alphay 0.5 \
    --output beckmann_walter07_0.5_0.5.png

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib ggx \
    --alphax 0.1 \
    --alphay 0.1 \
    --output ggx_walter07_0.1_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib ggx \
    --alphax 0.5 \
    --alphay 0.1 \
    --output ggx_walter07_0.5_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis false \
    --distrib ggx \
    --alphax 0.5 \
    --alphay 0.5 \
    --output ggx_walter07_0.5_0.5.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib beckmann \
    --alphax 0.1 \
    --alphay 0.1 \
    --output beckmann_heitz14_0.1_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib beckmann \
    --alphax 0.5 \
    --alphay 0.1 \
    --output beckmann_heitz14_0.5_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib beckmann \
    --alphax 0.5 \
    --alphay 0.5 \
    --output beckmann_heitz14_0.5_0.5.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib ggx \
    --alphax 0.1 \
    --alphay 0.1 \
    --output ggx_heitz14_0.1_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib ggx \
    --alphax 0.5 \
    --alphay 0.1 \
    --output ggx_heitz14_0.5_0.1.png

./build/microfacet \
    --samples 500 \
    --samplevis true \
    --distrib ggx \
    --alphax 0.5 \
    --alphay 0.5 \
    --output ggx_heitz14_0.5_0.5.png
