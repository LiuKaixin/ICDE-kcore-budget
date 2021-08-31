# ICDE-kcore-budget
Code Contributor: Kaixin Liu

If you have any questions, feel free to contact me. My email is lkx17@mails.tsinghua.edu.cn.

Please cite our paper if you choose to use our code. 

```
@inproceedings{DBLP:conf/icde/Liu0ZX21,
  author    = {Kaixin Liu and
               Sibo Wang and
               Yong Zhang and
               Chunxiao Xing},
  title     = {An Efficient Algorithm for the Anchored k-Core Budget Minimization
               Problem},
  booktitle = {37th {IEEE} International Conference on Data Engineering, {ICDE} 2021,
               Chania, Greece, April 19-22, 2021},
  pages     = {1356--1367},
  publisher = {{IEEE}},
  year      = {2021},
  url       = {https://doi.org/10.1109/ICDE51399.2021.00121},
  doi       = {10.1109/ICDE51399.2021.00121},
  timestamp = {Tue, 29 Jun 2021 12:06:35 +0200},
  biburl    = {https://dblp.org/rec/conf/icde/Liu0ZX21.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

## Tested Environment
- Ubuntu
- C++ 11
- GCC 4.8
- Boost
- cmake

## Compile
```sh
$ cmake .
$ make
```

## Parameters
```sh
./ICDE_kcore --operation query --algo <algorithm> [options]
```

- algo: which algorithm you prefer to run
    - exact: The exact algorithm
    - lkx: our algorithm
- options
    - --prefix \<prefix\>
    - --dataset \<dataset\>
    - -k \<degree constraint k\>
    - -b \<budget\>
    - -t \<threshold\>
    - --result_dir \<directory to place results\>

## Data
You can download from https://snap.stanford.edu/data/.



- Example

```sh
$ ./ICDE_kcore --algo clock --prefix ./data/ --dataset loc-gowalla -b 0.5 -k 10 -t 0.2
```



