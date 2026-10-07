[33mcommit 84ce8df64f74dffcba7437f041c109000b95b4bb[m[33m ([m[1;36mHEAD[m[33m -> [m[1;32mchunked[m[33m, [m[1;31morigin/chunked[m[33m, [m[1;32mmpi-chunked-backup[m[33m)[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Aug 27 12:14:10 2026 +0100

    putative mpi changes..

[33mcommit fa86918f88817879cdb271e4af3ce396769bc05c[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Jul 3 18:09:52 2026 +0100

    start mpi implementation

[33mcommit fc7150d47fed150f1e6f1751dbbcb104915419f2[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Jun 30 22:07:02 2026 +0100

    improve chunks

[33mcommit 36f67db78ba2d76de2326592bf66476f4e0e46e6[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Jun 30 20:48:03 2026 +0100

    improve status

[33mcommit 4bddd178e0381291e91c612ec6d98a1a0b3e2b99[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Jun 30 20:32:12 2026 +0100

    better parallel stuff

[33mcommit 2c3b4717bb6029c9c12c65eb7370919da0e2cda5[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Jun 29 20:57:29 2026 +0100

    clean up dot product implementation

[33mcommit 07294a79685cebe295f343a214ede4bb1c4a8aa8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Jun 29 20:34:42 2026 +0100

    working fast inner product

[33mcommit 3c19258af3f28e994f411b2b7aa47de91b27552f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Jun 26 12:34:05 2026 +0100

    automatically choose chunk size

[33mcommit 6d18224643f33db156bc8e9ed816489ad51bbaa9[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Jun 26 12:28:22 2026 +0100

    calculate chunks more efficiently

[33mcommit 6d189d1d35fe5ea6c17a1d02dedf646aa834769f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Jun 25 21:44:09 2026 +0100

    combinations

[33mcommit f8db8fcb9cf9a4161752b4ceeea5b94187e927f2[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Jun 25 11:27:16 2026 +0100

    improve output

[33mcommit 77d149a62245a3d471534f83b199e46652cbaaed[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Jun 24 13:11:35 2026 +0100

    attempt omp

[33mcommit 9701aa7ddcfa1121f4b61a21d1e8938343f8da06[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Jun 24 12:21:28 2026 +0100

    introduce chunks

[33mcommit 79e2c08b50293ef146785f864e1b862352843012[m[33m ([m[1;31morigin/main[m[33m, [m[1;31morigin/HEAD[m[33m, [m[1;32mmain[m[33m)[m
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Sun Jun 14 17:38:51 2026 +0100

    Batched gaussian elimination (#6)
    
    * Anticipate batched gaussian elimination
    
    * first effort
    
    * half of it done
    
    * working now
    
    * remove allocation
    
    * comment on reduced_error and potential fix for set_row
    
    * use combine instead of copy row
    
    * rref for a
    
    * row reduction on A
    
    * Do A row ops in place
    
    * Do row ops in place without needing A
    
    * Delete some superfluous code
    
    * remove some unnecessary parameters
    
    * reformat
    
    * use m4ri
    
    * revert dot product and reverse order of multiplication of errors
    
    * add option for fault counter
    
    * Fix multiplication order and s / k confusion
    
    * multiply destructively on left and right
    
    * changes
    
    * time to fault counter
    
    * option for fault counter cutoff
    
    * change multiplication order and add dry run

[33mcommit 46b2f59b3d9f2fb5fccd1ac3b6f2f7da66f84eec[m
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Tue May 12 10:39:28 2026 +0100

    Workstation fixes (#5)
    
    * some fixes
    
    * more fixes
    
    * revert some unnecessary changes

[33mcommit 9e69ef3e3304c2aecc36c38b8319d43ddc971419[m
Author: Callum Macpherson <93673602+CalMacCQ@users.noreply.github.com>
Date:   Tue Apr 28 12:07:32 2026 +0100

    Set up uv Python project, add basic output error script. (#3)
    
    * Set up a uv project, basic script showing usage of pauli-hugr
    
    * readme placeholder
    
    * explanatory comment
    
    * two -> 2 for consistency
    
    * update comments some more
    
    * allow Python 3.11, move files into src
    
    * update readme
    
    * use ssh and recommend uv sync --verbose
    
    * use tket 0.13

[33mcommit 42999834698b9824f6a9bbcacee51c6f28c15934[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 26 19:00:06 2026 +0100

    default log level

[33mcommit 605e0626404eb69e9d3d6ef005a8b46b03e67c68[m
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Sun Apr 26 18:59:27 2026 +0100

    Update README.org

[33mcommit daa4d1c0062e9c68889e84a9599827b66a61cd04[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 26 18:34:59 2026 +0100

    describe how to configure TIMING_MODE and LOG_LEVEL at build time

[33mcommit 0a79f981a644a912e10ceb997d3b65eb3c85f7a1[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 26 16:27:10 2026 +0100

    use clang-format

[33mcommit 3487bbde1f1e61ccd3b7b82d43d9ea3a0b2c74d6[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 26 15:56:54 2026 +0100

    format

[33mcommit aeb395f1542d7e255d83a0781a22810562206984[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 22:54:36 2026 +0100

    pretty up code some more

[33mcommit c2b7a278ddfd8f510612b7f26f606ec1040387e5[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 22:51:22 2026 +0100

    pretty up code

[33mcommit ba016030486a55504bf233c9bd8375bce170dcb2[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 22:44:07 2026 +0100

    remove some comments

[33mcommit 4b28c400a63887916093f05958906842c552e62f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 22:38:54 2026 +0100

    remove some comments

[33mcommit 70874f5da9e12c713486b5ae1abbd1beb57c399f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 20:12:38 2026 +0100

    clean up header file

[33mcommit 767ced61695f0f6ddb4ad59a1bdcad43da39815e[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 09:29:24 2026 +0100

    add thread count parameter

[33mcommit d142a8e9e1fdbaf4ce23c6e773382b364b4ea26f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 09:24:09 2026 +0100

    reindent

[33mcommit 1d24187ff060b797b72a64eba51d903d9365c070[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 00:32:25 2026 +0100

    fix coset verification

[33mcommit 2080ce32b326eeb42569b603df6ad36648097c2d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 25 00:23:10 2026 +0100

    update readme with list of potential improvements

[33mcommit 133994963c9fd0e6826f3d5673d6d864605fc566[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 23:04:22 2026 +0100

    improve build proces

[33mcommit 8b67cef81278580c5f780579fe38217766217287[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 22:25:24 2026 +0100

    clean up build process

[33mcommit 68c41d7054974b3e7b0f28216cc530b6d734bbd0[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 22:18:06 2026 +0100

    remove more superfluous code

[33mcommit dc8a95bc73e7733afcb8bca80fb8a137d85435c8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 22:12:29 2026 +0100

    remove superfluous code

[33mcommit 5ec1450ffb65205d408b598e67b9f8ebceba6440[m
Merge: df271c5 c688698
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 21:59:24 2026 +0100

    Merge pull request #2 from Quantinuum/correct-gaussian-elimination
    
    Correct gaussian elimination

[33mcommit c688698a6a321721336c2136cd3d32a564a3691c[m[33m ([m[1;31morigin/correct-gaussian-elimination[m[33m, [m[1;32mcorrect-gaussian-elimination[m[33m)[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 21:55:18 2026 +0100

    minor cleanups

[33mcommit 003c02701714fd6ceac3e07bbe1b54d0f6b1df78[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 17:19:21 2026 +0100

    maybe works now

[33mcommit df271c517497d98af0f1ad252e620c5cea0ff262[m
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 15:33:30 2026 +0100

    changes

[33mcommit 97618431a0e88487fd741e95f660b44b81587528[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 10:47:26 2026 +0100

    make

[33mcommit baa44f1111464914195809b82eedb9b6b2a66a13[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 10:08:28 2026 +0100

    another fix

[33mcommit dd3e0be14193fe456877b4ce3e1a9cbd26dfd1de[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 10:06:19 2026 +0100

    nearly working

[33mcommit 501f698c76db4d96784770639fc835e25d2b5a6e[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 24 09:51:16 2026 +0100

    start to fix gaussian elimination

[33mcommit 8061f6b6bbf19ee0688ad9eb866164146c1de163[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 15:53:34 2026 +0100

    pretty up the code

[33mcommit 13a2f06523b9be2dd85142a81f2e16360287b902[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 15:44:48 2026 +0100

    clean up c++

[33mcommit 28ceaac9850b4e1c6e2ffec6752675c134af0854[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 14:59:57 2026 +0100

    free memory when done with it

[33mcommit 116db30c9155bb74d437301e7cc6f17e93f2af2a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 14:52:56 2026 +0100

    allocate bit vectors a bit better

[33mcommit 47f1ac5cae4c6150205f2d66fceb7592bced8ea3[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 14:40:06 2026 +0100

    use calloc instead of malloc, safer

[33mcommit 9a4c3ecc0ed1b65ed7469402ae77a282e79de118[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 14:32:35 2026 +0100

    example usage from python

[33mcommit 3c976d1aa5710eb626be311edbe6eba1e327233a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 14:02:12 2026 +0100

    clean up python ffi

[33mcommit 1d73c67a228d23fe80bf3a124c54f986c1248711[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:56:54 2026 +0100

    turn off timing mode

[33mcommit be30bcef86a9f3a5e1b4d58a02a45f5fd76df8a1[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:50:30 2026 +0100

    update build instructions

[33mcommit ff0087b563f18ee8a949f6c7acda6bd78a33f9bf[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:48:49 2026 +0100

    use default c compiler and link to licence

[33mcommit 9f19296ad7ebd121bc1c5bf525e9755a6f0a6c6b[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:41:27 2026 +0100

    Improve python FFI

[33mcommit 3f9aa2aea77a0c8bc2410fc7c133fa8664216fa0[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:29:57 2026 +0100

    start FFI

[33mcommit 157cab31bbf238b90181484fadbc05f30f6e7191[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 23 13:09:54 2026 +0100

    reorganise project

[33mcommit c1947adf5d5b01902387755ffe4cf6e0d38a1f7e[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 17:01:31 2026 +0100

    fix testing

[33mcommit 637369142e46b0c0464428e3150bcb3865dfa17d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 16:02:09 2026 +0100

    another change

[33mcommit f19ab3f1962556e716694bb60d606c9130fe8a17[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 15:58:53 2026 +0100

    Update README.org

[33mcommit 2a107942b09921abf7a71a263fcb9eb63c8ebf3f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 15:56:00 2026 +0100

    fix makefile

[33mcommit 78df5cbf1a176bde914a3c28c6e460b41f7037fb[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 15:49:16 2026 +0100

    reorganise project tree

[33mcommit c414a81db8b13ea73b4a9033c176e35c2c148f38[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 22 15:21:04 2026 +0100

    build m4ri on demand

[33mcommit 4d4cd1d526cd7b48a79a4dd60bc29a1d00222098[m
Merge: f62945f 5104988
Author: selwynsimsek <selwyn.simsek@cantab.net>
Date:   Tue Apr 21 18:07:47 2026 +0100

    Merge pull request #1 from Quantinuum/dcp
    
    overhaul

[33mcommit 5104988530f990ec10c7bdbe9fc7ee663156e7ec[m[33m ([m[1;31morigin/dcp[m[33m, [m[1;32mdcp[m[33m)[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Apr 21 18:03:51 2026 +0100

    update readme

[33mcommit 1e18e1bee89f368aeb7c658dfeaa69d0253b345c[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Apr 21 11:44:51 2026 +0100

    reorganise tests

[33mcommit 69d4824dfeeb7820374f44793934abc66498f01d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 19 21:24:58 2026 +0100

    remove some superfluous code and add timing information

[33mcommit 8047e41ef280087a6342fa125f448c56ff4dd451[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 19 20:51:48 2026 +0100

    reorg

[33mcommit ce99c57423db8b399370aac57fbc7ac7af396c18[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 23:05:14 2026 +0100

    use popcount

[33mcommit a321587bb13a2a192251b18cbf9520a289bfa304[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 22:54:11 2026 +0100

    less transposing

[33mcommit 834734cb680f3c25528badb1c02467f3d7342f64[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 22:44:47 2026 +0100

    passes test suite??? think there are still some problems

[33mcommit 863afcc56bd758fc5276a2c6db0b534bc12af50a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 22:22:09 2026 +0100

    fixes

[33mcommit d0db24b71eaa136a9e8d3da028b3ad2c19fa6b1b[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 22:15:59 2026 +0100

    changes

[33mcommit 7437a3341210157d6af1fea01cf620e78be9bf54[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 21:28:30 2026 +0100

    gaussian elimination change

[33mcommit c568c350ebaf3decb62dc33519ee134575259370[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 21:17:41 2026 +0100

    improve debug

[33mcommit 1fe8203e36229a55a18ff6bb99a16d0468d08f7d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 21:07:05 2026 +0100

    new gaussian elimination

[33mcommit 1e3ed5dafcb5eeec1f418525355e965626b77fcf[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 20:49:34 2026 +0100

    change gaussian elimination

[33mcommit e2215a6fdc8dc22639f6bc35a22bd1faeb35f3fa[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 19:38:04 2026 +0100

    some changes

[33mcommit 39b19fae789fe15f9dee5b0e5bd02a37ccddddb1[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 17:17:20 2026 +0100

    clean up memory

[33mcommit 206f027cf5c952a53ef9360fe073dde93ca800d5[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 15:05:52 2026 +0100

    bit permutations

[33mcommit 099845f91473338ceeec595a7314c9da0dece261[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 14:25:45 2026 +0100

    rewrote dcp_fault_bad_p_css

[33mcommit 2138ff29dbd53c0c5507ecfaf976a1dfe4b266f0[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sat Apr 18 13:20:05 2026 +0100

    gaussian elimination / row ops

[33mcommit 0680a504fbbddd72c0380555fc3d324ea3d64bf8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Fri Apr 17 15:30:00 2026 +0100

    tests

[33mcommit f933b3effcd291823bca41e5392e6008381f308f[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 16 22:18:59 2026 +0100

    printing and other changes

[33mcommit 931f6ea62c75898e57650ddbbc4559918cf2c599[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 16 17:38:02 2026 +0100

    more changes

[33mcommit bbc8bbcce3ac8faebb5322782b6f6b85565b9610[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 16 14:21:30 2026 +0100

    moving to new conventions

[33mcommit 410801266ef1b05d5a53f3c02b4d0cf93a689e55[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Apr 14 12:21:33 2026 +0100

    try to ignore docs

[33mcommit aaa3e27cac33d81a1aeb2eb280255e1d1a9702e9[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Apr 14 12:20:24 2026 +0100

    more docs

[33mcommit 0a67f8a6f5445fc95a0b0ae295d24a92757d5c81[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Apr 14 11:50:48 2026 +0100

    start documenting

[33mcommit 9be8b520a7d8ff440598bcb9c5ebbec7c5d95a0d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 19:42:35 2026 +0100

    clean up code base

[33mcommit 8234a384a7c9a00cc4fab828aec7241a68ff4e82[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 13:07:42 2026 +0100

    another fix

[33mcommit e6352a2e0d18028819ac83c5046cb3fc36262dbd[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 13:03:51 2026 +0100

    fix

[33mcommit 8ebdcb81540e2f39a8e6eec480011639c36fa4a0[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 12:48:19 2026 +0100

    changes

[33mcommit 5aa0184b2dd5872a345c7889ad1fd89d6acf355c[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 12:08:54 2026 +0100

    remove permutation

[33mcommit 8b5d52d85a03d0d26d10d1f59241ac2bc95f6ed8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Apr 13 12:01:53 2026 +0100

    it works?

[33mcommit 3f7fe85a1c273739c22baabe226e7d4e4e88cf9d[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 21:57:28 2026 +0100

    changes

[33mcommit 7527b88fc3e25f969cca1947e0077632debe8cd9[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 20:32:47 2026 +0100

    some changes

[33mcommit 9b52e192ae08b435a00950d307d50e1251724769[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 20:00:05 2026 +0100

    some changes

[33mcommit a88313378400ab8fa012f10caffe925e511f8837[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 18:53:35 2026 +0100

    doing printing

[33mcommit 70c17c6dcb93006abaa83ac7753d391172cee35a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 09:57:31 2026 +0100

    more changes

[33mcommit e9b4c4002e0a9ee43b2fba2ec4f814d04aad3883[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 12 09:31:23 2026 +0100

    rref and permutation

[33mcommit 4a6f9015357bd44a77bd15dcbd489ded0e74a208[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 22:05:11 2026 +0100

    somehow working

[33mcommit fc894af47d3cf11c396d0aa9f18577f2820fffce[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 19:49:59 2026 +0100

    somehow done

[33mcommit 0dc96ad1ef51b234aa5024d02945783fd92150f1[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 18:29:25 2026 +0100

    start to do ft checking properly

[33mcommit b9e5a1d839df21de989670444b3e87ebf2858cae[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 16:43:40 2026 +0100

    finish converting

[33mcommit b3e41a2ed1096b3bcab9b7d0a5836df257e3ae75[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 14:31:58 2026 +0100

    more changes

[33mcommit 28048c0687fa68f2e601da56800abe4e8e2c625c[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 9 13:17:32 2026 +0100

    changes to api

[33mcommit ec258c33f782f6abd0c56029c3e8d10ca262c230[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 5 21:09:43 2026 +0100

    conversion nearly working

[33mcommit d5ed0c8a7ed6f8020bf2ee4c6083eeee22237741[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 5 19:13:17 2026 +0100

    changes

[33mcommit 3cf051822c1be875e04f6778066ab4f29591d1cd[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 5 16:25:14 2026 +0100

    last change before reshuffle

[33mcommit c5618530d705340ab2e8afeb8f86a6db19bddc6e[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Sun Apr 5 16:01:00 2026 +0100

    new interface

[33mcommit c9062dc989242534ff97623e71470196b18ea0d5[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Apr 2 11:50:01 2026 +0100

    commit before pushing

[33mcommit bd09c2d2164bf2874f9530ba97ef336600fba9bd[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 1 17:52:52 2026 +0100

    fix

[33mcommit a41157e0b7d39c8cf6f39cbbb05bcedd3854f646[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Apr 1 17:51:55 2026 +0100

    lots of changes

[33mcommit 9e513361954f43905bb46e78ca5e7f5376f361cc[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Mar 30 14:03:48 2026 +0100

    change dcp

[33mcommit 2ea84a5146b222d77c38c9a9063ec38497080709[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Mar 26 14:41:54 2026 +0000

    some more changes

[33mcommit d36a7c1e4eecfc5bc609a18f48d3097041272b9a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Mar 24 17:35:04 2026 +0000

    tentatively works

[33mcommit 78a4c69a563e73e2d64f05bded595a71b24f084a[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Mar 24 17:27:27 2026 +0000

    fix error in clifford multiplication

[33mcommit e9d88c6e26d326961fd7dc72301c6d6759687c76[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Mar 24 16:59:31 2026 +0000

    somewhat working but not really

[33mcommit b0e6d9bc12a57045ca8b730e98d6bd3c51f112a8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Tue Mar 24 15:36:34 2026 +0000

    changes

[33mcommit d3b81cbc1cf1e33c281bc6c200ce7f99363565c3[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Mar 18 16:31:49 2026 +0000

    work towards implementing diagonal clifford structure theorem

[33mcommit d55d818a6764eb39f34298b490afde5711712c31[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Mar 18 14:43:26 2026 +0000

    refactoring, start to implement diagonal clifford structure theorem

[33mcommit ae5c814159987541549db74090e52d672774eec8[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Mon Mar 16 23:36:44 2026 +0000

    some more routines

[33mcommit 92777f932dd1c6262f91555dadf3200508bdbec4[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Jan 22 10:50:58 2026 +0000

    init and free

[33mcommit e03c97b06d716c43908986efd51cb8ae1aa2b592[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Jan 22 10:14:10 2026 +0000

    new structs for paulis, cliffords and tableaux

[33mcommit bc53bf0c5c7419e1254d149e2e41727722e9afc3[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Jan 21 17:16:04 2026 +0000

    remove build script

[33mcommit 3c1d0536615ec1d83200da4c846e4d37b003c232[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Jan 21 17:12:56 2026 +0000

    new makefile

[33mcommit e0aed4d84cce10b6ea267135f68c4be73f15a06e[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Thu Jan 15 15:00:10 2026 +0000

    ignore m4ri

[33mcommit 0e3136895eecc80c574ca8d7dfe2222c693d8aad[m
Author: selwyn <selwyn.simsek@cantab.net>
Date:   Wed Jan 14 13:11:27 2026 +0000

    New interface

[33mcommit f62945f6a32f032f5c7fa1289e0a6a3d83d917fa[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Mon Mar 24 14:44:00 2025 +0000

    Undefault argument in C file.

[33mcommit dcabf8d6bb48e04fa20b6c59163059033a986a98[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Mon Mar 24 14:14:51 2025 +0000

    More tests.

[33mcommit f91b1863b45f26990c553be6a7624a4dfe263ba1[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Mon Mar 24 13:54:08 2025 +0000

    Update README.

[33mcommit 80fce6416eade608b8d9d616ed2d73fba1145cb0[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Sun Mar 23 13:30:54 2025 +0000

    Implement main function.

[33mcommit fde4d39173e88e0f3a81ce0fe3acca19476fa7b2[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Sun Mar 23 13:03:31 2025 +0000

    Handle zero-size buckets.

[33mcommit 9b00c8079de4402e0def3c652661f288803db0ea[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Sat Mar 22 18:30:05 2025 +0000

    Add class to iterate over product.

[33mcommit 17719a31feae5674e058a9a25ce9e5c35ebcefe0[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Sat Mar 22 12:44:04 2025 +0000

    Add class to iterate over combinations.

[33mcommit 85ceb09dea33973645b2b8938b4e73f35d076c4a[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Fri Mar 21 15:18:33 2025 +0000

    Add argument for number of worker threads.

[33mcommit fae380fae41d25608df9204801f694fd01439e40[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Mon Mar 24 13:00:23 2025 +0000

    fixup

[33mcommit 645d7c7e2d8d6dc43de228dd325f1647b9fa9236[m
Author: Alec Edgington <alec.edgington@quantinuum.com>
Date:   Fri Mar 21 15:07:42 2025 +0000

    Initial commit with skeleton.
