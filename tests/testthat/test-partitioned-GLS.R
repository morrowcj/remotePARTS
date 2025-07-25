data(ndvi_AK10000)
df = ndvi_AK10000[seq_len(1000), ] # first 1000 rows

## create partition matrix
# dput output sample_partitions(nrow(df), npart = 3)
pm = structure(c(40L, 929L, 931L, 400L, 982L, 177L, 739L, 983L, 76L,
            990L, 88L, 274L, 422L, 552L, 393L, 769L, 745L, 97L, 365L, 631L,
            710L, 743L, 119L, 102L, 856L, 371L, 460L, 960L, 998L, 152L, 9L,
            969L, 174L, 66L, 133L, 589L, 585L, 503L, 728L, 374L, 576L, 367L,
            247L, 159L, 706L, 711L, 50L, 811L, 373L, 724L, 634L, 242L, 500L,
            505L, 498L, 278L, 855L, 566L, 591L, 172L, 995L, 243L, 746L, 824L,
            354L, 874L, 2L, 90L, 875L, 885L, 953L, 627L, 535L, 175L, 455L,
            713L, 240L, 665L, 662L, 210L, 810L, 970L, 860L, 997L, 668L, 430L,
            883L, 789L, 835L, 441L, 999L, 854L, 404L, 821L, 218L, 197L, 691L,
            894L, 768L, 531L, 804L, 682L, 320L, 583L, 730L, 546L, 909L, 484L,
            194L, 486L, 216L, 652L, 128L, 313L, 408L, 312L, 618L, 644L, 336L,
            916L, 21L, 737L, 612L, 886L, 428L, 308L, 333L, 714L, 800L, 518L,
            120L, 513L, 897L, 890L, 78L, 280L, 857L, 17L, 16L, 459L, 406L,
            968L, 527L, 472L, 101L, 252L, 432L, 625L, 680L, 246L, 733L, 852L,
            369L, 699L, 6L, 75L, 742L, 735L, 501L, 499L, 734L, 81L, 809L,
            493L, 596L, 443L, 822L, 643L, 490L, 793L, 763L, 178L, 254L, 847L,
            707L, 433L, 851L, 966L, 812L, 181L, 586L, 613L, 352L, 55L, 626L,
            759L, 115L, 273L, 257L, 563L, 255L, 427L, 207L, 758L, 79L, 685L,
            941L, 603L, 623L, 621L, 315L, 823L, 414L, 449L, 951L, 831L, 321L,
            139L, 965L, 985L, 287L, 413L, 755L, 392L, 624L, 46L, 629L, 761L,
            145L, 358L, 979L, 777L, 87L, 363L, 187L, 778L, 642L, 420L, 219L,
            645L, 818L, 179L, 562L, 402L, 993L, 165L, 919L, 889L, 958L, 399L,
            157L, 186L, 471L, 111L, 984L, 438L, 697L, 8L, 142L, 672L, 359L,
            858L, 749L, 575L, 375L, 236L, 765L, 395L, 454L, 339L, 752L, 540L,
            158L, 271L, 826L, 632L, 656L, 22L, 692L, 13L, 223L, 224L, 108L,
            125L, 548L, 925L, 550L, 865L, 23L, 84L, 107L, 660L, 987L, 412L,
            927L, 185L, 636L, 696L, 382L, 285L, 893L, 519L, 14L, 558L, 649L,
            86L, 895L, 981L, 217L, 694L, 263L, 901L, 173L, 774L, 262L, 124L,
            48L, 581L, 725L, 27L, 25L, 237L, 434L, 57L, 418L, 42L, 190L,
            944L, 437L, 361L, 282L, 322L, 568L, 106L, 994L, 494L, 954L, 192L,
            819L, 516L, 720L, 648L, 123L, 671L, 11L, 744L, 536L, 837L, 135L,
            928L, 950L, 283L, 600L, 248L, 164L, 738L, 512L, 667L, 602L, 451L,
            940L, 838L, 411L, 551L, 146L, 547L, 485L, 917L, 15L, 286L, 18L,
            1L, 182L, 196L, 450L, 524L, 232L, 149L, 504L, 291L, 674L, 98L,
            249L, 56L, 905L, 806L, 241L, 446L, 844L, 85L, 892L, 934L, 391L,
            195L, 52L, 880L, 760L, 416L, 340L, 162L, 136L, 325L, 543L, 109L,
            292L, 574L, 815L, 307L, 708L, 41L, 841L, 775L, 372L, 62L, 904L,
            445L, 896L, 53L, 229L, 368L, 492L, 978L, 873L, 848L, 617L, 276L,
            495L, 537L, 639L, 731L, 155L, 461L, 869L, 683L, 659L, 71L, 961L,
            814L, 33L, 370L, 132L, 764L, 532L, 99L, 377L, 469L, 19L, 988L,
            522L, 921L, 211L, 463L, 570L, 794L, 675L, 390L, 959L, 640L, 996L,
            456L, 915L, 483L, 816L, 957L, 304L, 409L, 560L, 226L, 166L, 718L,
            442L, 356L, 986L, 279L, 256L, 297L, 510L, 171L, 29L, 309L, 717L,
            314L, 289L, 284L, 92L, 323L, 153L, 732L, 381L, 638L, 405L, 366L,
            876L, 299L, 261L, 497L, 955L, 318L, 561L, 1000L, 137L, 805L,
            573L, 657L, 930L, 541L, 341L, 328L, 756L, 67L, 104L, 834L, 467L,
            779L, 222L, 347L, 220L, 529L, 239L, 212L, 462L, 112L, 715L, 64L,
            598L, 474L, 332L, 300L, 750L, 238L, 594L, 903L, 389L, 72L, 344L,
            555L, 35L, 59L, 798L, 888L, 593L, 787L, 666L, 771L, 693L, 188L,
            571L, 130L, 266L, 569L, 839L, 797L, 557L, 465L, 938L, 845L, 863L,
            193L, 183L, 180L, 943L, 799L, 935L, 864L, 712L, 539L, 10L, 397L,
            68L, 268L, 549L, 39L, 853L, 786L, 740L, 884L, 388L, 825L, 654L,
            980L, 877L, 597L, 80L, 305L, 690L, 65L, 337L, 508L, 82L, 349L,
            70L, 189L, 453L, 114L, 69L, 355L, 937L, 447L, 661L, 203L, 167L,
            380L, 134L, 528L, 898L, 678L, 491L, 275L, 45L, 294L, 350L, 346L,
            870L, 7L, 425L, 741L, 475L, 424L, 669L, 956L, 701L, 118L, 487L,
            73L, 103L, 947L, 913L, 306L, 900L, 622L, 419L, 792L, 417L, 709L,
            515L, 926L, 448L, 327L, 559L, 121L, 949L, 747L, 784L, 542L, 705L,
            60L, 54L, 605L, 871L, 204L, 489L, 20L, 199L, 517L, 716L, 906L,
            545L, 727L, 801L, 387L, 473L, 77L, 96L, 650L, 767L, 817L, 436L,
            464L, 403L, 310L, 722L, 228L, 590L, 49L, 963L, 881L, 592L, 479L,
            974L, 936L, 633L, 989L, 398L, 607L, 421L, 971L, 781L, 595L, 113L,
            100L, 588L, 407L, 206L, 991L, 36L, 878L, 466L, 977L, 673L, 523L,
            26L, 565L, 511L, 131L, 599L, 850L, 105L, 201L, 933L, 298L, 319L,
            862L, 213L, 200L, 700L, 245L, 882L, 584L, 117L, 807L, 702L, 353L,
            401L, 202L, 975L, 676L, 250L, 110L, 601L, 684L, 480L, 386L, 233L,
            783L, 468L, 227L, 564L, 866L, 939L, 721L, 664L, 234L, 899L, 28L,
            317L, 221L, 533L, 376L, 338L, 335L, 647L, 138L, 184L, 163L, 431L,
            820L, 260L, 141L, 751L, 635L, 231L, 766L, 842L, 423L, 230L, 689L,
            351L, 44L, 496L, 394L, 288L, 911L, 74L, 579L, 520L, 867L, 530L,
            37L, 31L, 3L, 215L, 198L, 482L, 384L, 670L, 704L, 846L, 736L,
            426L, 615L, 577L, 168L, 521L, 160L, 776L, 829L, 478L, 364L, 908L,
            770L, 914L, 832L, 973L, 830L, 458L, 38L, 723L, 967L, 942L, 782L,
            235L, 470L, 58L, 147L, 122L, 43L, 687L, 270L, 918L, 4L, 924L,
            302L, 32L, 879L, 34L, 628L, 932L, 677L, 127L, 440L, 572L, 488L,
            383L, 922L, 126L, 553L, 457L, 251L, 514L, 507L, 608L, 868L, 902L,
            641L, 95L, 148L, 208L, 225L, 891L, 964L, 176L, 833L, 646L, 83L,
            976L, 357L, 754L, 803L, 920L, 620L, 385L, 609L, 962L, 653L, 360L,
            51L, 477L, 326L, 827L, 796L, 757L, 91L, 506L, 910L, 439L, 269L,
            378L, 303L, 861L, 209L, 567L, 972L, 330L, 679L, 259L, 611L, 663L,
            859L, 265L, 89L, 214L, 582L, 410L, 795L, 30L, 610L, 61L, 316L,
            604L, 150L, 296L, 415L, 244L, 923L, 481L, 170L, 5L, 791L, 556L,
            526L, 295L, 616L, 534L, 293L, 587L, 872L, 681L, 651L, 331L, 753L,
            329L, 544L, 788L, 140L, 525L, 205L, 301L, 538L, 698L, 253L, 554L,
            719L, 156L, 396L, 780L, 695L, 637L, 12L, 324L, 808L, 785L, 348L,
            688L, 116L, 509L, 703L, 435L, 948L, 93L, 946L, 630L, 272L, 129L,
            992L, 773L, 476L, 762L, 726L, 345L, 343L, 342L, 772L, 840L, 144L,
            802L, 258L, 290L, 887L, 729L, 614L, 658L, 619L, 191L, 151L, 47L,
            334L, 606L, 311L, 267L, 945L, 452L, 264L, 849L, 379L, 836L, 502L,
            362L, 580L, 154L, 578L, 912L, 429L, 169L, 655L, 277L, 686L, 143L,
            813L, 63L, 161L, 843L, 748L, 281L, 907L, 952L, 790L, 444L, 94L,
            828L), dim = c(333L, 3L), dimnames = list(NULL, c("part.1", "part.2",
                                                              "part.3")))

## fit GLS with fixed nugget
partGLS = fitGLS_partition(formula = CLS_coef ~ 0 + land, partmat = pm,
                           data = df, nugget = 0, do.t.test = TRUE, progressbar = FALSE)

## hypothesis tests
chisqr(partGLS) # explanatory power of model
t.test(partGLS) # significance of predictors

## now with a numeric predictor
partGLS_2 <- fitGLS_partition(formula = CLS_coef ~ lat, partmat = pm, data = df, nugget = 0,
                              progressbar = FALSE)

test_that("sample_partitions works properly", {
  # 20 random samples from 1 to 1e5
  tmp = sample_partitions(npix = 1e5, npart = 20)
  expect_equal(dim(tmp), c(1e5/20, 20))
  expect_true(all(tmp %in% 1:1e5))
  # random samples of size 3000 from 1 to 1e5
  tmp2 = sample_partitions(npix = 1e5, npart = NA, partsize = 3000)
  expect_equal(dim(tmp2), c(3000, 1e5 %/% 3000))
  expect_true(all(tmp2 %in% 1:1e5))
  # random samples of size 5000 from 1e3 to 1e5
  tmp3 = sample_partitions(pixels = 1e3:1e5, partsize = 5000, npart = NA)
  expect_equal(dim(tmp3), c(5000, length(1e3:1e5) %/% 5000))
  expect_true(all(tmp3 %in% 1e3:1e5))
})

test_that("partitioned GLS $part works correctly for categorical predictors", {
  expect_equal(partGLS$part$coefficients,
               structure(c(-0.000143365351062626, 9.04337643353381e-05, 2.96698091401277e-06,
                           -0.000441376628359838, -0.00114798051178247, -0.00127815627178324,
                           -0.000146605758337422, -3.00018157420377e-05, 0.000352513590903927),
                         dim = c(3L, 3L), dimnames = list(NULL, c("landShrubland",
                                                           "landSavanna", "landGrassland"))),
               tolerance = 1e-3)
  expect_equal(partGLS$part$covar_coefs,
               structure(c(3.15102787130498e-06, 3.02342945641165e-06, 2.98185691282179e-06,
                           3.02342945641164e-06, 3.19312250487829e-06, 2.95756275475129e-06,
                           2.98185691282179e-06, 2.95756275475129e-06, 3.08575788170884e-06,
                           3.85251577802123e-06, 3.77386836022111e-06, 3.68242877039011e-06,
                           3.77386836022111e-06, 3.96755406841202e-06, 3.65660140926314e-06,
                           3.68242877039011e-06, 3.65660140926314e-06, 3.75855418673139e-06,
                           2.25109645111239e-06, 2.12801934327663e-06, 2.06157559818322e-06,
                           2.12801934327663e-06, 2.25329802916544e-06, 2.03761893639001e-06,
                           2.06157559818322e-06, 2.03761893639001e-06, 2.25840489416393e-06
               ), dim = c(3L, 3L, 3L),
               dimnames = list(c("landShrubland", "landSavanna", "landGrassland"),
                               c("landShrubland", "landSavanna", "landGrassland"), NULL)),
               tolerance = 1e-8)
  expect_equal(partGLS$part$pvals_t,
               structure(c(0.935678558089016, 0.963278928321389, 0.998423374784677,
                           0.805059795291598, 0.564783348646881, 0.39512002432872, 0.933537648722641,
                           0.987662396868634, 0.814686923835494), dim = c(3L, 3L),
                         dimnames = list(NULL, c("landShrubland", "landSavanna", "landGrassland"))),
               tolerance = 1e-3)
})

test_that("partitioned GLS $part works correctly for numerical predictors", {
  expect_equal(partGLS_2$part$coefficients,
               structure(c(-0.0219162336196283, -0.0056201970601145, -0.0314005971712678,
                           0.000344487039671403, 8.75587717251568e-05, 0.000486740443703192
               ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "lat"))),
               tolerance = 1e-5)
  expect_equal(partGLS_2$part$covar_coefs,
               structure(c(0.000394897213286253, -6.22178670917822e-06, -6.22178670917804e-06,
                           9.87809325991826e-08, 0.000376483705374476, -6.00810980389132e-06,
                           -6.00810980389129e-06, 9.68528029189814e-08, 0.000446868840384066,
                           -6.9560703556472e-06, -6.95607035564721e-06, 1.08808011498555e-07
               ), dim = c(2L, 2L, 3L), dimnames = list(c("(Intercept)", "lat"
               ), c("(Intercept)", "lat"), NULL)),
               tolerance = 1e-6)
  expect_equal(partGLS_2$part$pvals_t,
               structure(c(0.2708853062044, 0.772262854830753, 0.138385625352877,
                           0.273847134737324, 0.778619430284168, 0.141003412152979), dim = 3:2, dimnames = list(
                             NULL, c("(Intercept)", "lat"))),
               tolerance = 1e-3)
})

test_that("partitioned GLS $cross works correctly for categorical predictors", {
  expect_equal(partGLS$cross$rcoefs,
               structure(c(0.965498600896112, 1.03954390770187, 1.10539039997512,
                           0.955010802332451, 1.02888984612658, 1.08885451266607, 0.965951196095157,
                           1.02741337189767, 1.09547277476428, 0.950329067869466, 1.03298811392531,
                           1.0988481772292, 0.950604399967651, 1.03345984048255, 1.09509493107161,
                           0.953509881257318, 1.02519466514889, 1.09256023644712, 0.962564040073074,
                           1.02845109262656, 1.08958141410839, 0.95446614649029, 1.02019308961727,
                           1.07554707771412, 0.974006113134597, 1.0309820522337, 1.09621242959237
               ), dim = c(3L, 3L, 3L), dimnames = list(NULL, c("landShrubland", "landSavanna", "landGrassland"),
                                                       c("landShrubland", "landSavanna", "landGrassland"))),
               tolerance = 1e-4)
  expect_equal(partGLS$cross$rSSRs,
               c(0.0202893512228277, 0.0141380004492578, 0.0240324517112229),
               tolerance = 1e-4)
  expect_equal(partGLS$cross$rSSEs,
               c(0.27964680867385, 0.25297944660295, 0.235174561342431),
               tolerance = 1e-4)
})

test_that("partitioned GLS $cross works correctly for numerical predictors", {
  expect_equal(partGLS_2$cross$rcoefs,
               structure(c(0.853434239080802, 0.846797232191893, 0.793344390171755,
                           -0.862979733845368, -0.846104047033159, -0.790114764722982, -0.851342138627944,
                           -0.856806525933185, -0.811690916793361, 0.869648072580667, 0.86272247072584,
                           0.816529730938086), dim = c(3L, 2L, 2L),
                         dimnames = list(NULL, c("(Intercept)", "lat"), c("(Intercept)", "lat"))),
               tolerance = 1e-4)
  expect_equal(partGLS_2$cross$rSSRs, c(0.756287770143775, 0.744290061495832, 0.666720801505638),
               tolerance = 1e-4)
  expect_equal(partGLS_2$cross$rSSEs, c(0.277899981881308, 0.250741531749591, 0.232404205436463),
               tolerance = 1e-4)
})

test_that("partitioned GLS $overall works correctly for categorical predictors", {
  expect_equal(partGLS$overall$coefficients,
               c(landShrubland = -1.66548686044251e-05, landSavanna = -0.000955837803975185,
                 landGrassland = 5.86353389414892e-05),
               tolerance = 1e-6)
  expect_equal(partGLS$overall$rcoefficients,
               structure(c(1.03681096952437, 1.02425172037503, 1.0296124475857,
                           1.02738845300799, 1.02638639050727, 1.02375492761777, 1.02686551560268,
                           1.01673543794056, 1.03373353165355), dim = c(3L, 3L), dimnames = list(
                             c("landShrubland", "landSavanna", "landGrassland"),
                             c("landShrubland", "landSavanna", "landGrassland"))),
               tolerance = 1e-4)
  expect_equal(partGLS$overall$t.test,
               structure(c(-1.66548686044251e-05, -0.000955837803975185, 5.86353389414892e-05,
                           0.00176725721179705, 0.00177527922842883, 0.00175198282994246,
                           -0.00942413390266458, -0.538415472151459, 0.0334679872081937,
                           0.992482633303081, 0.590410912296134, 0.973308119902519), dim = 3:4, dimnames = list(
                             c("landShrubland", "landSavanna", "landGrassland"),
                             c("Est", "SE", "t.stat", "pval.t"))),
               tolerance = 1e-5)
  expect_equal(partGLS$overall$Fstat, 2.46317202739637, tolerance = 1e-4)
})

test_that("partitioned GLS $overall works correctly for numerical predictors", {
  expect_equal(partGLS_2$overall$coefficients,
               c(`(Intercept)` = -0.0196456759503369, lat = 0.000306262085033251),
               tolerance = 1e-4)
  expect_equal(partGLS_2$overall$rcoefficients,
               structure(c(0.831191953814816, -0.83306618186717, -0.839946527118163,
                           0.849633424748197), dim = c(2L, 2L),
                         dimnames = list(c("(Intercept)", "lat"), c("(Intercept)", "lat"))),
               tolerance = 1e-4)
  expect_equal(partGLS_2$overall$t.test,
               structure(c(-0.0196456759503369, 0.000306262085033251, 0.0189719846603902,
                           0.000302078004379305, -1.03550979520626, 1.01385099409188, 0.300682238792957,
                           0.310900079692848), dim = c(2L, 4L),
                         dimnames = list(c("(Intercept)", "lat"), c("Est", "SE", "t.stat", "pval.t"))),
               tolerance = 1e-5)
  expect_equal(partGLS_2$overall$Fstat, 1.15263133122574, tolerance = 1e-4)
})
