Plot.Period = function(y, t.selected) {
    plot(y[t.selected] ~ t.selected, type = "l");
};

Periodic.Linear.Model = function(y) {
    t = 1:length(y);
    periodic.model = lm(y ~ t + I(sin(t * 2 * pi / 7)) + I(cos(t * 2 * pi / 7)));
    plot(y ~ t, type = "l");
    lines(fitted(periodic.model) ~ t, col = "blue");
    return(periodic.model);
};

Max.Claimable.Beta = function(periodic.linear.fit, error.level = 0.05, precision = 1e-4) {
    beta.OLS = coefficients(periodic.linear.fit)[2];
    beta.OLS.SE = sqrt(vcov(periodic.linear.fit)[2, 2]);
    beta.candidate = seq(from = 0, to = beta.OLS, by = precision);
    t.crit = qt(error.level, df = periodic.linear.fit$df.residual, lower.tail = F);
    for (beta in beta.candidate) {
        t.stat = (beta.OLS - beta) / beta.OLS.SE;
        if (t.stat <= t.crit) {
            return(beta - precision);
        }
    }
};
        
        

covid.df = read.csv("owid-covid-data.csv", header = T);
covid.SA.df = covid.df[covid.df$location == "South Africa",];
covid.SA.date = covid.SA.df$date;
covid.SA.date[1];
covid.SA.date[length(covid.SA.date)];
I.new = covid.SA.df$new_cases;
# The data for South Africa has a number of NAs at the start, and some daily readings
# are zero. We remove the preceding NAs and replace intermediate zeros and NAs with
# interpolated values.
I.new[1:50];
covid.SA.date = covid.SA.date[41:(length(covid.SA.date) - 1)];
I.new = I.new[41:(length(I.new) - 1)];
t = 1:length(covid.SA.date);
covid.SA.date[1];
covid.SA.date[length(covid.SA.date)];
length(I.new);
plot(I.new ~ t, type = "l");
ln.I.new = log(I.new);
sum(is.infinite(ln.I.new));
for (i in 1:length(ln.I.new)) {
    if (is.infinite(ln.I.new[i]) | is.na(ln.I.new[i])) {
        ln.I.new[i] = (ln.I.new[i - 1] + ln.I.new[i + 1]) / 2;
    }
}
plot(ln.I.new ~ t, type = "l", ylab = "ln(R(t))", xlab = "t",
     main = "Trend of ln(R(t))) for South Africa's covid transmission");

t1 = 30:100;
t2 = 250:290;
t3 = 410:470;
t.omicron = 614:length(t);
length(t.omicron);
covid.SA.date[t1];
covid.SA.date[t2];
covid.SA.date[t3];
covid.SA.date[t.omicron];

# Estimated beta for episode 1.
plot(ln.I.new[t1], type = "l");
t1.fit = Periodic.Linear.Model(ln.I.new[t1]);
summary(t1.fit);
shapiro.test(residuals(t1.fit));
beta1 = t1.fit$coefficients[2];
SE.beta1 = sqrt(vcov(t1.fit)[2, 2]);

# Estimated beta for episode 2.
Plot.Period(ln.I.new, t2);
t2.fit = Periodic.Linear.Model(ln.I.new[t2]);
summary(t2.fit);
shapiro.test(residuals(t2.fit));
beta2 = coefficients(t2.fit)[2];
SE.beta2 = sqrt(vcov(t2.fit)[2, 2]);

# Estimated beta for episode 3.
Plot.Period(ln.I.new, t3);
t3.fit = Periodic.Linear.Model(ln.I.new[t3]);
summary(t3.fit);
shapiro.test(residuals(t3.fit));
beta3 = coefficients(t3.fit)[2];
SE.beta3 = sqrt(vcov(t3.fit)[2, 2]);

# Estimated beta for the omicron episode.
Plot.Period(ln.I.new, t.omicron);
t.omicron.fit = Periodic.Linear.Model(ln.I.new[t.omicron]);
summary(t.omicron.fit);
shapiro.test(residuals(t.omicron.fit));
beta.omicron = coefficients(t.omicron.fit)[2];
SE.beta.omicron = sqrt(vcov(t.omicron.fit)[2, 2]);

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2));
plot(ln.I.new[t1] ~ t1, ylab = "residual", xlab = "t", main = "Episode I");
lines(fitted(t1.fit) ~ t1, col = "blue");
plot(ln.I.new[t2] ~ t2, ylab = "residual", xlab = "t", main = "Episode II");
lines(fitted(t2.fit) ~ t2, col = "blue");
plot(ln.I.new[t3] ~ t3, ylab = "residual", xlab = "t", main = "Episode III");
lines(fitted(t3.fit) ~ t3, col = "blue");
plot(ln.I.new[t.omicron] ~ t.omicron, ylab = "residual", xlab = "t",
     main = "Omicron Episode");
lines(fitted(t.omicron.fit) ~ t.omicron, col = "blue");

start = covid.SA.date[c(t1[1], t2[1], t3[1], t.omicron[1])];
end = c(covid.SA.date[c(t1[length(t1)], t2[length(t2)], t3[length(t3)])], NA);
beta = c(beta1, beta2, beta3, beta.omicron);
SE.beta = c(SE.beta1, SE.beta2, SE.beta3, SE.beta.omicron);
CI1 = beta1 + c(-1, 1) * qt(0.025, df = t1.fit$df.residual, lower.tail = F) * SE.beta1;
CI2 = beta2 + c(-1, 1) * qt(0.025, df = t2.fit$df.residual, lower.tail = F) * SE.beta2;
CI3 = beta3 + c(-1, 1) * qt(0.025, df = t3.fit$df.residual, lower.tail = F) * SE.beta3;
CI.omicron = beta.omicron + c(-1, 1) * qt(0.025, df = t.omicron.fit$df.residual, lower.tail = F) *
    SE.beta.omicron;
parameter.df = data.frame(
    Episode = c("I", "II", "III", "Omicron"),
    Start = start,
    End = end,
    Beta = round(beta, 4),
    SE.Beta = round(SE.beta, 5),
    Beta.Max = c(Max.Claimable.Beta(t1.fit), Max.Claimable.Beta(t2.fit),
                 Max.Claimable.Beta(t3.fit), Max.Claimable.Beta(t.omicron.fit)),
    CI.Left = round(c(CI1[1], CI2[1], CI3[1], CI.omicron[1]), 4),
    CI.Right = round(c(CI1[2], CI2[2], CI3[2], CI.omicron[2]), 4)
);
parameter.df;

k.candidate = seq(from = 1, to = beta.omicron / beta3, by = 0.01);
z.crit = qnorm(0.05, mean = 0, sd = 1, lower.tail = F);
for (k in k.candidate) {
    z.stat = (beta.omicron - k * beta3) / sqrt((SE.beta.omicron^2) + (SE.beta3^2) * (k^2));
    if (z.stat <= z.crit) {
        print(k - 0.01);
        break;
    }
}

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2));
plot(residuals(t1.fit) ~ t1, ylab = "residual", xlab = "t", main = "Episode I");
plot(residuals(t2.fit) ~ t2, ylab = "residual", xlab = "t", main = "Episode II");
plot(residuals(t3.fit) ~ t3, ylab = "residual", xlab = "t", main = "Episode III");
plot(residuals(t.omicron.fit) ~ t.omicron, ylab = "residual", xlab = "t",
     main = "Omicron Episode");
