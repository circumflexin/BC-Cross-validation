


# UAV19_253_1
require(boot)
BAI1 = c(9.60399449,
9.701547531,
9.802350427,
9.736363636,
9.733382031)

n = 10000
boot_mean <- function(original_vector, resample_vector) {
  mean(original_vector[resample_vector])
}
BD_boots1 = boot(BAI1, boot_mean, n)
plot(density(BD_boots1$t))

#UAV19_254_2
BAI2 = c(BAI1,
9.935944363)
median(BAI2)

n = 10000
BD_boots2 = boot(BAI2, boot_mean, n)
plot(density(BD_boots2$t))
plot(density(BAI2))
boot.ci(BD_boots2,type="bca")
median(BAI2)


# UAV19_251_1
BAI3 = c(9.573845153,
9.197869593,
10.2388852,
10.17027864,
9.716061828,
9.558020478)

n = 10000

BD_boots3 = boot(BAI3, boot_mean, n)
plot(density(BD_boots3$t))
plot(density(BAI3))
median(BAI3)
boot.ci(BD_boots3,type="bca")

