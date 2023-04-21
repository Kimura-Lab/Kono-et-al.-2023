# AHP (Analytic Hierachy Process) による評価を行う
AHP <- function(	x,					# 評価基準の重要度（下三角行列をベクトルで用意）
			y,					# 代替案の評価（各行列の下三角行列を列とする行列で用意）
			labels.x=NULL,				# 評価基準のラベル
			labels.y=NULL)				# 代替案のラベル
{
	items <- function(n)					# 下三角行列の要素数から行列サイズを求める
	{
		retval <- (1+sqrt(1+8*n))/2
		return(if (retval!=floor(retval)) Inf else retval)
	}
	make.matrix <- function(x)				# 正方行列から重みベクトルを求める
	{
		n <- items(length(x))				# 行列のサイズ
		mat <- diag(n)					# 三角行列を表すベクトルから行列を生成
		mat[lower.tri(mat, diag=FALSE)] <- x
		mat <- t(mat)+mat
		mat[upper.tri(mat)] <- 1/mat[upper.tri(mat)]
		diag(mat) <- 1
		result <- eigen(mat)				# 固有値・固有ベクトルを求める
		val <- as.numeric(result$values[1])
		vec <- as.numeric(result$vectors[,1])
		weight <- vec/sum(vec)				# 固有ベクトルを和が 1 になるように標準化したものが重み
		ci <- (val-n)/(n-1)
		cr <- ci/c(0,0,0.58,0.9,1.12,1.24,1.32,1.41,1.45,1.49,1.51,1.53)[n]
		if (ci > 0.1 || cr > 0.1) {
			cat("\nC.I.=", ci, ",  C.R.=", cr, "\n", sep="")
			print(mat)
			W <- outer(weight, weight, "/")
			print(W)
			print(mat-W)
		}
		return(list(lambda=val, vec=vec, weight=weight, ci=ci, cr=cr))
	}
	if (is.null(labels.x)) {				# ラベルが与えられていないときは A, B, ...
		labels.x <- LETTERS[1:items(length(x))]
	}
	ans.x <- make.matrix(x)
	weight.x <- ans.x$weight				# 評価基準の重要度
	names(weight.x) <- labels.x
	nitems.y <- items(nrow(y))
	if (is.null(labels.y)) {				# ラベルが与えられていないときは a, b, ...
		labels.y <- letters[1:nitems.y]
	}
	ans.y <- matrix(unlist(apply(y, 2, make.matrix)), 3+2*nitems.y)
	weight.y <- ans.y[(2+nitems.y):(1+2*nitems.y),]		# 代替案の評価
	rownames(weight.y) <- labels.y
	colnames(weight.y) <- labels.x
	score <- rowSums(t(weight.x*t(weight.y)))		# スコア
	return(structure(list(weight.x=weight.x, weight.y=weight.y, score=score, sorted.score=sort(score)), class="AHP"))
}
# print メソッド
print.AHP <- function(	obj,					# AHP の返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	cat("\n評価基準の重み\n\n")
	print(round(obj$weight.x, digits=digits))
	cat("\n代替案の評価結果\n\n")
	print(round(obj$weight.y, digits=digits))
	cat("\nスコア\n\n")
	print(round(obj$score, digits=digits))
	cat("\nソートされたスコア\n\n")
	print(round(obj$sorted.score, digits=digits))
}
# plot メソッド
plot.AHP <- function(	obj,					# AHP の返すオブジェクト
			xlab="Score",				# 結果グラフの横軸名
			main="AHP (Analytic Hierachy Process)",	# 結果グラフの表題
			file="")				# 結果グラフをファイル出力するときにファイル名
{
	if (file != "") pdf(file, width=540/72, height=160/72, onefile=FALSE)
	score <- obj$score
	plot(score, rep(0, length(score)), pch=19, xlab=xlab, main=main, xaxt="n",
		xlim=range(pretty(score)), ylab="", yaxt="n", ylim=c(0,0.2),
		bty="n", xpd=TRUE)
	text(score, 0.0, names(score), pos=3)
	axis(1, pos=0)
	if (file != "") dev.off()
}
# AIC による，ヒストグラム（度数分布表）の最適階級分割の探索
#
AIC.Histogram <- function(	x,					# データベクトル
				d = 0,					# 測定精度（無限の精度の場合には 0
				c = floor(2*sqrt(length(x))-1)) {	# 初期階級数
	MODEL <- function(c1, r, c2, n, c)				# 個別の度数分布表の AIC の計算
	{
		logNZ <- function(x, y) {				# 補助関数
			return(ifelse(x > 0, x*log(x/y) , 0))
		}
		N <- sum(n)						# サンプルサイズ
		sum1 <- sum(n[1:c1])					# 左端で併合され階級について
		sum2 <- sum(n[(c-c2+1):c])				# 右端で併合され階級について
		cc1c2r1 <- (c-c1-c2)/r+1				# 中央部分で併合され階級について
		temp <- 0
		for (j in 2:cc1c2r1) {
			sum3 <- sum(n[(c1+(j-2)*r+1):(c1+(j-1)*r)])
			temp <- temp+logNZ(sum3, r*N)
		}
		AIC <- -2*(logNZ(sum1, c1*N)+temp+logNZ(sum2, c2*N))+2*cc1c2r1 # モデルの AIC
		return(AIC)
	}
	x <- x[!is.na(x)]						# 完全なデータについて
	N <- length(x)							# サンプルサイズ
	brks = seq(min(x)-d/2, max(x)+d/2, length=c+1)			# 分割点
	ans <- hist(x, breaks=brks, plot=FALSE, right=FALSE)		# 度数分布表を得る
	n <- ans$count							# 度数ベクトル
	result <- NULL							# 結果の蓄積用
	for (c1 in 1:(c-2)) {						# 左端で併合する階級数を探索
		for (c2 in 1:(c-2)) {					# 右端で併合する階級数を探索
			if (c <= c1+c2) next				# 制約条件を満たさない場合は次へ
			for (r in 1:(c-c1-c2)) {			# 中央付近で併合する階級数を探索
				if ((c-c1-c2)%%r == 0) {		# 併合する階級数は等しくする必要がある
					result <- append(result, list(c1, r, c2, MODEL(c1, r, c2, n, c)))
				}
			}
		}
	}
	p <- n/N*100							# パーセント
	df <- data.frame(matrix(unlist(result), ncol=4, byrow=TRUE))	# 結果をデータフレームに変換
	colnames(df) <- c("c1", "r", "c2", "AIC")			# 列に名前を付ける
	o <- order(df[,4])						# AIC の小さい順
	df <- df[o,]							# 並べ替える
	plot(range(brks), c(0, max(p)), type="n", xlab="data", ylab="Percent")		# プロット枠組み
	sapply(1:(N-1), function(i) rect(brks[i], 0, brks[i+1], p[i], border="gray"))	# 初期ヒストグラム
	c1 <- df[1,1]							# 最適モデルの結果取り出し
	r <- df[1,2]
	c2 <- df[1,3]
	AIC <- df[1,4]
	lo <- brks[1]							# 階級の開始値
	delta <- diff(brks[1:2])					# 階級幅
	m <- (c-c1-c2)/r						# 中央付近の併合後階級数
	cmg <- cumsum(c(0, c1, rep(r, m), c2))				# 併合後の階級の開始値ベクトル
	sapply(1:(m+2), function(i) rect(lo+cmg[i]*delta, 0,            # 併合後のヒストグラム
		lo+cmg[i+1]*delta, mean(p[(cmg[i]+1):(cmg[i+1])]),
		border="red", col="red", density=15))
	mtext(paste("AIC =", round(AIC, 2)), side=3, line=-1.2,         # AIC を図に書き込む
		at=lo+0.8*c*delta)
	invisible(list(df=df, n=n, breaks=brks))			# 結果を返す
}
# AIC による分割表の独立性の判定
AIC.independence <- function(x)
{
	log2 <- function(n) sum(ifelse(n == 0, 0, n*log(n)))	# Σ n・log(n)

	rt <- rowSums(x)					# 行和
	ct <- colSums(x)					# 列和
	n <- sum(x)						# 総和（合計）
	lr <- nrow(x)						# 行数
	lc <- ncol(x)						# 列数
	AIC0 <- 2*(2*log2(n)-log2(rt)-log2(ct)+lr+lc-2)		# 二要因が独立であるとするモデルの AIC
	AIC1 <- 2*(log2(n)-log2(x)+lr*lc-1)			# 二要因が独立でないとするモデルの AIC
	result <- ifelse(AIC0 < AIC1, "二要因は独立である", "二要因は独立ではない") 
	list(AIC.independent=AIC0, AIC.dependent=AIC1, result=result)
}
# 二要因の分散分析（ASB タイプ；SPFp・q デザイン；混合計画）
ASB <- function(data)								# 3次元配列のデータ
{										# 次元は，被験者，要因B，要因A の順
	n <- dim(data)[1]							# 水準の組み合わせでの被験者数
	nm1 <- n-1
	Nb <- dim(data)[2]							# 要因 B の水準の数
	Na <- dim(data)[3]							# 要因 A の水準の数
	grand.mean <- mean(data)						# 全体の平均
	e.a <- sum((apply(data, 3, mean)-grand.mean)^2)*Nb*n			# 要因A
	diff.obj <- sum((apply(data, c(1, 3), mean)-grand.mean)^2)*Nb-e.a	# 個人差（S）
	e.b <- sum((apply(data, 2, mean)-grand.mean)^2)*Na*n			# 要因B
	cross <- sum((apply(data, c(3, 2), mean)-grand.mean)^2)*n-e.a-e.b	# 交互作用 A x S
	err <- sum(apply(data, c(3, 2), var)*nm1)-diff.obj			# 交互作用 S x B
	SS <- c(e.a, diff.obj, e.b, cross, err)					# 平方和
	dfa <- Na-1
	dfb <- Nb-1
	df <- c(dfa, Na*nm1, dfb, dfa*dfb,  Na*nm1*dfb)				# 自由度
	MS <- SS/df								# 平均平方
	f <- p <- rep(NA, 5)
	f[c(1, 3, 4)] <- MS[c(1, 3, 4)]/MS[c(2, 5, 5)]				# F 値
	p[c(1, 3, 4)] <- pf(f[c(1, 3, 4)], df[c(1, 3, 4)], df[c(2, 5, 5)],	# P 値
			    lower.tail=FALSE)
	result <- data.frame(SS, df, MS, f, p)
	colnames(result) <- c("SS", "d.f.", "MS", "F value", "P value")
	rownames(result) <- c("Factor A", "S", "Factor B", "AxS", "SxB")
	class(result) <- c("anova.table", "data.frame")
	return(result)
}
# 総当たり法による重回帰分析を行う
# 　データフレームには，分析に使用する独立変数と従属変数のみを含むこと。
# 　また，従属変数は最終列に置くこと。
#
All.possible.subset.selection <- function(df,				# データフレーム（独立変数，従属変数）
					  limit=10)			# 独立変数の個数の上限（数が多いと計算時間が指数的に増える）
{
	df <- subset(df, complete.cases(df))				# 欠損値を持つケースを除く
	nv <- ncol(df)-1						# 独立変数の個数
	if (nv > limit) {						# limit より多いと分析を中断する
		stop(paste("独立変数が", limit,
		"個以上である（多すぎる）。\n",
		"limit 引数で変更できる", paste=""))
	}
	n <- 2^nv							# 独立変数を取り出す取り出し方
	bincomb <- matrix(FALSE, nrow=n, ncol=nv)			# e1071 パッケージの bincombinations より
	for (j in 1:nv) {
		bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
	}
	bincomb <- bincomb[-1,]
	n <- n-1
	name <- names(df)						# 変数名を取り出す
	depname <- name[nv+1]
	name <- name[1:nv]
	result4 <- character(n)						# 数値型ベクトル確保
	result1 <- result2 <- result3 <- numeric(n)			# 数値型ベクトル確保
	for (i in 1:n) {						# 独立変数の全ての組み合わせについて，
		str <- name[bincomb[i,]]				# どの独立変数が使われるかを割り出す
		form <- reformulate(str, depname)			# モデル式を作る（"formula" クラス）
		ans <- lm(form, df)					# 重回帰分析の結果
		result <- summary(ans)
		result1[i] <- result$r.square				# 重相関係数の二乗（決定係数）
		result2[i] <- result$adj.r.square			# 自由度調整済み重相関係数の二乗
		result3[i] <- AIC(ans)					# AIC
		temp <- as.character(form)				# モデル式を文字列に変換
		result4[i] <- paste(temp[2], "~", temp[3])		# モデル式を記録
	}
	return(structure(list(rsq=result1, adj=result2, AIC=result3, form=result4),
	class="all.possible.subset.selection"))
}
print.all.possible.subset.selection <- function( obj,			# "all.possible.subset.selection" クラスのオブジェクトをプリント
					 	 sort.by=c("adj", "rsq", "AIC"), # 結果を何で並べ替えるかを指示
						 models=20)	 	# 良い方から何番目まで出力するか
{
	result <- data.frame(obj$rsq, obj$adj, obj$AIC, obj$form)
	sort.by <- match.arg(sort.by)
	o <- order(switch (sort.by, "rsq"=result[,1], "adj"=result[,2], "AIC"=result[,3]), decreasing = sort.by != "AIC")
	result <- result[o,]
	cat("\nR square  Adjusted R square       AIC       Formula\n")	# 表頭
	models <- min(models, nrow(result))
	apply(result[1:models,], 1, function(x)				# 各行の出力
		cat(sprintf("%8.5f %13.5f  %13.5f       %s\n",
		as.double(x[1]), as.double(x[2]), as.double(x[3]), x[4])))
	invisible(result)
}
# 総当たり法による重回帰分析を行う
# 　データフレームには，分析に使用する独立変数と従属変数のみを含むこと。
# 　また，従属変数は最終列に置くこと。
#
All.possible.subset.selection2 <- function(  df,			# データフレーム（独立変数，従属変数）
						   fixed,		### 追加：必ず取り入れる変数の名前（ベクトル）
						   limit=10)		# 独立変数の個数の上限（数が多いと計算時間が指数的に増える）
{
	df <- subset(df, complete.cases(df))			# 欠損値を持つケースを除く
	nv <- ncol(df)-1-length(fixed)					### 変更　独立変数の個数
	if (nv > limit) {						# limit より多いと分析を中断する
		stop(paste("独立変数が", limit,
		"個以上である（多すぎる）。\n",
		"limit 引数で変更できる", paste=""))
	}
	n <- 2^nv							# 独立変数を取り出す取り出し方
	bincomb <- matrix(FALSE, nrow=n, ncol=nv)			# e1071 パッケージの bincombinations より
	for (j in 1:nv) {
		bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
	}
	bincomb <- bincomb[-1,]
	n <- n-1
	name <- names(df)						# 変数名を取り出す
	depname <- name[ncol(df)]					### 変更
	name <- setdiff(name[1:(ncol(df)-1)], fixed)		### 変更
	result4 <- character(n)					# 数値型ベクトル確保
	result1 <- result2 <- result3 <- numeric(n)		# 数値型ベクトル確保
	for (i in 1:n) {						# 独立変数の全ての組み合わせについて，
		str <- c(name[bincomb[i,]], fixed)			### 変更　どの独立変数が使われるかを割り出す
		form <- reformulate(str, depname)			# モデル式を作る（"formula" クラス）
		ans <- lm(form, df)					# 重回帰分析の結果
		result <- summary(ans)
		result1[i] <- result$r.square			# 重相関係数の二乗（決定係数）
		result2[i] <- result$adj.r.square			# 自由度調整済み重相関係数の二乗
		result3[i] <- AIC(ans)				# AIC
		temp <- as.character(form)				# モデル式を文字列に変換
		result4[i] <- paste(temp[2], "~", temp[3])	# モデル式を記録
	}
	return(structure(list(rsq=result1, adj=result2, AIC=result3, form=result4),
	class="all.possible.subset.selection"))
}
print.all.possible.subset.selection <- function( obj,		# "all.possible.subset.selection" クラスのオブジェクトをプリント
					 	 sort.by=c("adj", "rsq", "AIC"), # 結果を何で並べ替えるかを指示
						 models=20)	 	# 良い方から何番目まで出力するか
{
	result <- data.frame(obj$rsq, obj$adj, obj$AIC, obj$form)
	sort.by <- match.arg(sort.by)
	o <- order(switch (sort.by, "rsq"=result[,1], "adj"=result[,2], "AIC"=result[,3]), decreasing = sort.by != "AIC")
	result <- result[o,]
	cat("\nR square  Adjusted R square       AIC       Formula\n")	# 表頭
	models <- min(models, nrow(result))
	apply(result[1:models,], 1, function(x)			# 各行の出力
		cat(sprintf("%8.5f %13.5f  %13.5f       %s\n",
		as.double(x[1]), as.double(x[2]), as.double(x[3]), x[4])))
	invisible(result)
}
# 三角多項式グラフ（Andrews グラフ）
Andrews.graph <- function(	dat,				# データ行列
				normalize=TRUE,			# 変数ごとに正規化する
				points=100,			# 曲線のなめらかさを決める
				col=NULL,
				...)				# plot への引数
{
	dat <- subset(dat, complete.cases(dat))
	dat <- as.matrix(dat)					# 行列に変換
	if (normalize == TRUE) {				# normalize が TRUE なら，
		dat <- scale(dat)				# 変数ごとに正規化する
	}
	n <- nrow(dat)						# 行数
	nv <- ncol(dat)						# 列数（変数の数）
	t <- seq(-pi, pi, length=points)			# -π ～ πに points 個の点を取る
	coef <- matrix(sapply(1:nv,				# 係数行列
		function(i) if (i %% 2 == 0) sin((i%/%2)*t)
			    else cos((i%/%2)*t)),
		nrow=nv, byrow=TRUE)
	coef[1,] <- rep(1/sqrt(2), points)
	data <- matrix(sapply(1:n, function(k) colSums(dat[k,]*coef)), byrow=TRUE, nr=n)
	plot(range(t), range(data), type="n",
		xaxt="n", xlab="", ylab="")
	axis(1, at=seq(-pi, pi, length=5),
	     lab=expression(-pi, -pi/2, 0, pi/2, pi))
	if (is.null(col)) {
		col <- rep("blue", n)
	}
	else if (length(col) == 1) {
		col <- rep(col, n)
	}
	invisible(sapply(1:n, function(k) lines(t, data[k,], col=col[k])))
}
# Breslow-Day 検定
BD.test <- function(m)		# 2×2×k の配列
{
	method <- "Breslow-Day 検定"
	data.name <- deparse(substitute(m))
	k <- dim(m)[3]
	Nk <- apply(m, 3, sum)
	psiMH <- sum(m[1, 1,]*m[2, 2,]/Nk)/sum(m[2, 1,]*m[1, 2,]/Nk)
	nk1 <- m[1, 1,]+m[1, 2,]
	nk2 <- m[2, 1,]+m[2, 2,]
	xkt <- m[1, 1,]+m[2, 1,]
	a1 <- psiMH-1
	b1 <- -psiMH*(nk1+xkt)-nk2+xkt
	c1 <- psiMH*nk1*xkt
	e <- (-b1-sqrt(b1^2-4*a1*c1))/(2*a1)
	v <- 1/(1/e+1/(nk1-e)+1/(xkt-e)+1/(nk2-xkt+e))
	chisqBD <- sum((m[1, 1,]-e)^2/v)
	df <- k-1
	p <- pchisq(chisqBD, df, lower.tail=FALSE)
	return(structure(list(statistic=c("chi sq."=chisqBD), parameter=c(df=df), p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# ブラッドリー・テリーのモデル
# 東京大学教養学部統計学教室編「基礎統計学III　自然科学の統計学」東京大学出版会
# によったが，他のところにあるのとは答えが違う
BTM <- function(x,						# 一対比較の結果を表す正方行列
		constant=1,					# 解の制約条件（解の合計値）
		max.rotation=500,				# 収束計算を行う上限回数
		epsilon=1e-10)					# 収束判定値
{
	data.name <- deparse(substitute(x))
	method <- "ブラッドリー・テリーのモデル"
	nc <- ncol(x)						# 項目数
	stopifnot(nc == nrow(x))				# 正方行列でないと分析中止
	if (is.null(dimnames(x))) {				# 項目名がないときは補完する
		labels <- LETTERS[1:nc]
	}
	else if(is.null(colnames(x))) {
		labels <- rownames(x)
	}
	else {
		labels <- colnames(x)
	}
	diag(x) <- 0
	r <- x+t(x)
	yi. <- rowSums(x)
	theta <- rep(constant/nc, nc)
	names(theta) <- labels
	for (i in 1:max.rotation) {
		theta <- yi./colSums(r/outer(theta, theta, "+"))
		theta.new <- theta*(constant/sum(theta))
		converge <- all(abs(theta.new-theta) < epsilon)
		theta <- theta.new
		if (converge) {
			break
		}
	}
	if (!converge) {
		warning("収束しませんでした。max.rotation を大きくして再実行してみてください")
	}
	s <- !diag(nc)
	expected <- (r*theta/outer(theta, theta, "+"))
	w <- expected[s]
	chi2 <- sum((x[s]-w)^2/w)
	df <- choose(nc-1, 2)
	P <- pchisq(chi2, df, lower.tail=FALSE)
	names(chi2) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=chi2, parameter=df, 
        	p.value=P, method=method, data.name=data.name, observed=x, 
        	expected=expected, theta=theta), class=c("htest", "BTM")))
}
# summary メソッド
summary.BTM <- function(obj,						# BTM オブジェクト
			digits=5)					# 表示する桁数
{
	cat("\nスコア\n\n")
	print(round(obj$theta, digits=digits))
	cat("\nソートされたスコア\n\n")
	print(round(sort(obj$theta), digits=digits))
	cat("\n観察値\n\n")
	print(round(obj$observed, digits=digits))
	cat("\n期待値\n\n")
	print(round(obj$expected, digits=digits))
}
# plot メソッド
plot.BTM <- function(	obj,						# BTM オブジェクト
			xlab="Score",					# 結果グラフの横軸名
			main="Bradley-Terry model (Paired Comparison)",	# 結果グラフの表題
			file="")					# 結果グラフをファイル出力するときにファイル名
{
	theta <- obj$theta
	if (file != "") pdf(file, width=540/72, height=160/72, onefile=FALSE)
	plot(theta, rep(0, length(theta)), pch=19, xlab=xlab, main=main, xaxt="n",
		xlim=range(pretty(theta)), ylab="", yaxt="n", ylim=c(0,0.2),
		bty="n", xpd=TRUE)
	text(theta, 0.0, names(theta), pos=3)
	axis(1, pos=0)
	if (file != "") dev.off()
}
# バートレットの球面性検定
Bartlett.sphericity.test <- function(x)                 # データ行列
{
       method <- "Bartlett's test of sphericity"
       data.name <- deparse(substitute(x))
       x <- subset(x, complete.cases(x))                # 欠損値を持つケースを除く
       n <- nrow(x)
       p <- ncol(x)
       chisq <- (1-n+(2*p+5)/6)*log(det(cor(x)))
       df <- p*(p-1)/2
       p.value <- pchisq(chisq, df, lower.tail=FALSE)
       names(chisq) <- "X-squared"
       names(df) <- "df"
       return(structure(list(statistic=chisq, parameter=df, p.value=p.value,
              method=method, data.name=data.name), class="htest"))
}
# ボンフェローニの方法および関連する手法による多重比較
Bonferroni <- function(	data,						# 観察値ベクトル
			group,						# 群変数ベクトル
			method=c("Bonferroni", "Holm", "Shaffer", "Holland-Copenhaver"),
			# ボンフェローニ，ホルム，シェーファー，ホランド・コペンハーバー
			alpha=0.05)					# 有意水準
{
	OK <- complete.cases(data, group)				# データと群変数がともに欠損値ではないデータを抽出する
	data <- data[OK]
	group <- factor(group[OK])
	method <- match.arg(method)					# method に与えられた略語による指定を補完する
# ボンフェローニ法　Bonferroni
	n <- table(group)						# 各群のデータ数
	a <- length(n)							# 群の数
	k <- a*(a-1)/2							# 2 群の総当たり数
	m <- tapply(data, group, mean)					# 各群の平均値
	v <- tapply(data, group, var)					# 各群の不偏分散
	phi <- length(data)-a						# 自由度
	V <- sum((n-1)*tapply(data, group, var))/phi			# 誤差分散
	t <- combn(a, 2, function(ij) abs(diff(m[ij]))/sqrt(sum(V/n[ij])))	# abs(m[i]-m[j])/sqrt(V*(1/n[i]+1/n[j])))
	P <- pt(t, phi, lower.tail=FALSE)*2				# P 値を計算
	result1 <- cbind("n"=n, "mean"=m, "variance"=v)			# 各群の，標本サイズ，平均値，不偏分散
	rownames(result1) <- paste("Group", 1:a, sep="")		# 表側を作る
	result2 <- cbind("t value"=t, "P value"=P)			# 二群の比較結果の，t 値，P 値
	rownames(result2) <- combn(a, 2, paste, collapse=":")		# 表側を作る
	if (method == "Bonferroni") {					# ボンフェローニ法 Bonferroni のとき，
		P.threshold <- alpha/k					# P 値の限界値と，
		judge <- result2[,2] <= P.threshold			# 判定（有意のときには 1，有意でないときには 0 が表示される）
	}
	else {								# ボンフェローニ法以外のとき
		result2 <- result2[order(result2[,2]),]
		if (method != "Holm" && a > 9) {			# ホルムの方法以外で，群数が 9 より大きいときには対応できないのでホルムの方法にする
			warning("Too many groups. Use Holm.")
			method = "Holm"
		}
		if (method == "Holm") {					# ホルムの方法 Holm のとき
			P.threshold <- alpha/(k:1)			# P 値の限界値を設定
		}
		else {
			h <- c(3,1,1, 6,rep(3,3),2,1, 10,rep(6,4),4,4:1, 15,rep(10,5),rep(7,3),6,4,4:1,
			       21,rep(15,6),rep(11,4),10,9,7,7:1, 28,rep(21,7),rep(16,5),15,13,13:1,
			       36,rep(28,8),rep(22,6),21,rep(18,3),16,16,15,13,13:1)[(a*(a-1)*(a-2)/6):((a^3-a)/6-1)]
			P.threshold <- if (method == "Shaffer") alpha/h else 1-(1-alpha)^(1/h)	# P 値の限界値を設定
		}
		judge <- cumprod(result2[,2] <= P.threshold) != 0	# 判定（有意のときには 1，有意でないときには 0 が表示される）
		if (method == "Holm" || method == "Shaffer" || method == "Holland-Copenhaver") {
			min.pos <- which.min(judge)			# 途中で有意でないものが出てきたら，それ以降は全て結果を保留
			if (judge[min.pos] == 0) {			# 初めて出てくる 0 の位置
				judge[min.pos:k] <- NA			# 保留を表すという意味で NA を代入
			}
		}
	}
	result2 <- cbind(result2,					# 結果の追加
		    "P threshold"=P.threshold, "judge"=judge)
	return(list(method=method, alpha=alpha, n.of.tests=k, result2=result2,
		remark="judge=1 significant, judge=0 not significant, judge=NA suspend", result1=result1, phi=phi, V=V))
}
# 分散・共分散行列の同等性の検定
BoxM <- function(	x,						# データ行列（データフレーム）
			gvar)						# 群を表す変数
{
	method <- "分散・共分散行列の同等性の検定"
	data.name <- paste(deparse(substitute(x)), "~", deparse(substitute(gvar)))
	x <- as.data.frame(x)
	nv <- ncol(x)							# 変数の個数
	if (nv < 2) stop("分散共分散行列を計算する変数は2個以上必要です")
	gvar <- as.factor(gvar)
	ni <- table(gvar)						# 各群のサンプルサイズ
	n <- length(gvar)						# サンプルサイズ
	g <- length(ni)							# 群の数
	y <- split(x, gvar)						# 群ごとに分割したデータ行列
	Si <- lapply(y, var)						# 分散共分散行列
	log.det.Si <- sapply(Si, function(x) log(det(x)))		# 行列式の対数値
	S <- sapply(y, function(x) (nrow(x)-1)*var(x))			# 変動行列
	S <- matrix(rowSums(S), nv, nv)/(n-g)				# プールした変動行列
	M <- (n-g)*log(det(S))-sum((ni-1)*log.det.Si)			# Box の M 統計量
	f1 <- (g-1)*nv*(nv+1)/2						# 第 1 自由度
	rho <- 1-(2*nv^2+3*nv-1)/(6*(nv+1)*(g-1))*(sum(1/(ni-1))-1/(n-g))
	tau <- (nv-1)*(nv+2)/(6*(g-1))*(sum(1/(ni-1)^2)-1/(n-g)^2)
	f2 <- (f1+2)/abs(tau-(1-rho)^2)					# 第 2 自由度
	gamma <- (rho-f1/f2)/f1
	F <- M*gamma							# F 値
	p <- pf(F, f1, f2, lower.tail=FALSE)				# P 値
	return(structure(list(statistic=c(M=M, F=F),
		parameter=c(df1=f1, df2=f2), p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# Box-Cox 変換の，最適のλを図により求める
Box.Cox.transformation <- function(	x,			# データ
					l=-3, r=3,		# λの探索範囲
					delta=0.1)		# 探索ステップ
{
	x <- x[!is.na(x)]					# 有効データのみを対象とする
	Gm <- exp(mean(log(x)))					# 幾何平均を求める
	lambda0 <- seq(l, r, delta)				# 候補とするλの値のベクトル
	result <- numeric(length(lambda0))			# 変換値を保存するベクトル
	for (i in seq(along=lambda0)) {				# λの各候補ごとに，
		lambda <- lambda0[i]
		if (lambda == 0) {
			w <- Gm*log(x)
		}
		else {
			w <- (x^lambda-1)/(lambda*Gm^(lambda-1))
		}
		result[i] <- sd(w)				# 変換結果を保管
	}
	plot(lambda0, result, type="l")				# 図に描く
}
# Box-Cox 変換の，最適のλをシンプレックス法によって求める
Box.Cox.transformation2 <- function(	x,				# データ
					loop = 500,			# 収束計算の許容ループ数
					epsilon = 1e-15,		# 収束判定値
					Alpha = 2,			# シンプレックス法のα
					Beta = 0.5,			# シンプレックス法のβ
					Gamma = 2)			# シンプレックス法のγ
{
	Box.Cox.sub <- function(lambda)					# ボックス・コックス変換の計算値
	{
		if (lambda == 0) {
			w <- Gm*log(x)
		}
		else {
			w <- (x^lambda-1)/(lambda*Gm^(lambda-1))
		}
		sd(w)
	}

	x <- x[!is.na(x)]						# 有効データのみを対象とする
	Gm <- exp(mean(log(x)))						# 幾何平均を求める

	p1 <- -3							# 初期値
	p2 <- p1+0.1							# 初期値

	vec <- c(p1, p2)

	for (i in 1:loop) {						# 収束計算

		result <- c(Box.Cox.sub(vec[1]), Box.Cox.sub(vec[2]))

		h <- which.max(result)					# 大きい方
		s <- which.min(result)					# 小さい方

		ph <- vec[h]
		fh <- result[h]
		ps <- vec[s]
		fs <- result[s]

		p0 <- vec[s]

		pr <- (1+Alpha)*p0-Alpha*ph
		fr <- Box.Cox.sub(pr)

		if (fr > fh) {						# 探索値の状況により推測値を修正
			pc <- Beta*ph+(1-Beta)*p0
			vec[h] <- pc
		}
		else if (fs > fr) {
			pe <- Gamma*pr+(1-Gamma)*p0
			fe <- Box.Cox.sub(pe)
			if (fr > fe) vec[h] <- pe else vec[h] <- pr
		}
		else {
		 	vec[h] <- pr
		}
		if (abs((vec[1]-vec[2])/vec[1]) < epsilon) break	# 近似値の差が小さいなら収束と見なす
	}

	mean(vec)							# 2つのベクトルの平均値は，真値により近いだろう
}
Brown.Forsythe.test <- function(x, group) {  # x: データベクトル，group: 群変数ベクトル
  data.name <- paste(deparse(substitute(x)), "~", deparse(substitute(group)))
  OK <- complete.cases(x, group)				# 欠損値を持つケースを除く
  x <- x[OK]
  group <- as.factor(group[OK])
  group <- group[, drop=TRUE]
  d <- split(x, group)
  df1 <- length(d)-1
  ni <- sapply(d, length)
  mi <- sapply(d, mean)
  ui <- sapply(d, var)
  ci <- (1-ni/sum(ni))*ui
  F.BF <- sum(ni*(mi-mean(x))^2)/sum(ci)
  C <- ci/sum(ci)
  df2 <- 1/sum(C^2/(ni-1))
  p <- pf(F.BF, df1, df2, lower.tail=FALSE)
  method <- "Brown-Forsythe 検定"
  return(structure(list(statistic=c(F=F.BF), parameter=c("df1"=df1, "df2"=df2), "p.value"=p, method=method, data.name=data.name), class="htest"))
}
# 中央値の（差の）信頼区間を求める
CI4median <- function(	x,						# データベクトル
			y = NULL,					# 一標本の場合には y は省略する
			conf.level = 0.95,				# 信頼率
			method = c(	"independent.sample",		# 独立二標本
					"paired.sample",		# 対応のある標本
					"one.sample"))			# 一標本
{
	stopifnot(conf.level > 0, conf.level < 1)			# 信頼率は割合で指定する
	method <- match.arg(method)					# 引数の補完
	if (is.null(y)) {						# y が NULL なら，
		method <- "one.sample"					# 一標本
	}

	if (method == "independent.sample") {				# 独立二標本（中央値の差の信頼区間）
		x <- x[!is.na(x)]					# 欠損値を持つケースを除く
		y <- y[!is.na(y)]					# 欠損値を持つケースを除く
		n1 <- length(x)						# サンプルサイズ
		n2 <- length(y)						# サンプルサイズ
		vec <- sort(as.vector(outer(x, y, "-")))		# 全ての組み合わせで引き算し，小さい順に並べる
		k <- ifelse(n1 <= 100 && n2 <= 100,			# 分位点を計算
			qwilcox(0.5-conf.level/2, n1, n2),
			n1*n2/2-qnorm(0.5+conf.level/2)*sqrt(n1*n2*(n1+n2+1)/12))
		n <- n1*n2						# 統計量の値の総数
	}
	else {
		if (method == "paired.sample") {			# 対応のある標本（の中央値の差の信頼区間）
			stopifnot(is.null(y) == FALSE)			# y にもデータがあるはず
			OK <- complete.cases(x, y)			# 欠損値を持つケースを除く
			x <- x[OK]
			y <- y[OK]
			x <- y-x					# 一標本の中央値の信頼区間に帰結する
		}
		n1 <- length(x)						# サンプルサイズ
		mean <- outer(x, x, "+")/2				# あらゆる組み合わせで平均値を取る
		vec <- sort(mean[upper.tri(mean, diag=TRUE)])		# 小さい順に並べる
		k <- ifelse(n1 <= 300,					# 分位点を計算
			qsignrank(0.5-conf.level/2, n1),
			n1*(n1+1)/4-qnorm(0.5+conf.level/2)*sqrt(n1*(n1+1)*(2*n1+1)/24))
		n <- n1*(n1+1)/2					# 統計量の値の総数
	}
	result <- c("LCL"=vec[k], "UCL"=vec[n+1-k])			# 信頼限界値
	return(list(name = method, result = result))
}
# コクラン・アーミテージ検定
Cochran.Armitage <- function(	r.i,			# 反応数のベクトル
				n.i,			# ケース数のベクトル
				x.i=seq(along=r.i))	# 外的基準。省略されたときは1から始まる整数値が仮定される
{
	k <- length(r.i)				# 群の数
	stopifnot(length(n.i) == k, length(x.i) == k)

	p.i <- r.i/n.i					# 標本比率のベクトル
	r <- sum(r.i)					# 反応数の合計
	n <- sum(n.i)					# ケース数の合計
	p.m <- r/n					# プールした標本比率（標本比率の平均値）
	x.m <- sum(n.i*x.i)/n				# 外的基準の平均値
	xx <- x.i-x.m					# 平均偏差
	b <- sum(n.i*(p.i-p.m)*xx)/sum(n.i*xx^2)	# 傾き
	a <- p.m-b*x.m					# 切片
	xt <- b^2*sum(n.i*xx^2)/(p.m*(1-p.m))		# 傾き（トレンド）に対するカイ二乗値
	xh <- n^2*(sum(r.i^2/n.i)-r^2/n)/r/(n-r)	# 非一様性
	xq <- xh-xt					# 直線からの乖離に対するカイ二乗値
	chisq <- c(xt, xq, xh)				# カイ二乗値
	df <- c(1, k-2, k-1)				# 自由度
	P <- pchisq(chisq, df, lower.tail=FALSE)	# P 値
	res <- data.frame(chisq=chisq, df=df, P=P)
	colnames(res) <- c("カイ二乗値", "自由度", "P 値")
	rownames(res) <- c("トレンド", "直線からの乖離", "非一様性")
	return(res)
}
# コクランの Q 検定
Cochran.Q.test <- function(x)					# データ行列
{
	data.name <- deparse(substitute(x))
	method <- "コクランの Q 検定"
	x <- subset(x, complete.cases(x))			# 欠損値を持つケースを除く
	k <- ncol(x)						# 変数の個数（列数）
	g <- colSums(x)						# 列和
	l <- rowSums(x)						# 行和
	Q <- ((k-1)*(k*sum(g^2)-sum(g)^2))/(k*sum(l)-sum(l^2))	# 検定統計量
	df <- k-1						# 自由度
	p <- pchisq(Q, df, lower.tail=FALSE)			# P 値
	names(Q) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=Q, parameter=df, p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# 星座グラフを描く
Constellation.graph <- function(
		dat,						# 多変量データ
		w=NULL,						# 変数ごとの重みベクトル
		points=100,					# 半円をなめらかに描くための点の数
		...)						# plot, point に渡す引数
{
	dat <- t(as.matrix(subset(dat, complete.cases(dat))))	# データ行列は転置したものを使う
	if (is.null(w)) {					# 重みベクトルのデフォルトは，
		n <- nrow(dat)					# 変数の個数で，
		w <- rep(1/n, n)				# 1 を割った値（同じ値）
	}
	mx <- apply(dat, 1, max)				# 変数ごとの最大値
	mn <- apply(dat, 1, min)				# 変数ごとの最小値
	rg <- mx-mn						# 変数ごとの範囲
	dat <- pi*(dat-mn)/rg					# 円内に収まるように正規化
	y <- colSums(sin(dat)*w)				# x 座標
	x <- colSums(cos(dat)*w)				# y 座標
	plot(c(-1, 1), c(0, 1), type="n",			# 描画枠を決める
		xaxt="n", xlab=" ", yaxt="n", ylab="", asp=1, ...)
	points(x, y, ...)					# 点を打つ
	axis(1, at=-1:1, lab=-1:1, pos=0)			# 横軸
	axis(2, at=0:1, lab=0:1)				# 縦軸
	theta <- seq(0, pi, length=points)			# 半円を描く準備
	lines(c(0, cos(theta), 0, 0), c(0, sin(theta), 0, 1)) 	# 半円を描く
}
# Cox-Mantel 検定
Cox.Mantel <- function(	group,					# 群を識別するベクトル（1, 2 のいずれか）
			event,					# 死亡なら 1，生存なら 0 の値をとるベクトル
			time)					# 生存期間ベクトル
{
	method <- "Cox-Mantel 検定"
	data.name <- sprintf("time: %s, event: %s, group: %s",
		deparse(substitute(time)), deparse(substitute(event)), deparse(substitute(group)))
	OK <- complete.cases(group, event, time)		# 欠損値を持つケースを除く
	group <- group[OK]
	event <- event[OK]
	time <- time[OK]
	stopifnot(length(group) == length(event),		# 要素数が同じであること
		  length(group) == length(time))
	len <- length(group)					# 要素数
	tg <- table(time, group*10+event)			# 群とエンドポイントごとに生存時間をまとめる
	k <- nrow(tg)						# 行数
	nia <- table(group)[1]					# 第 1 群のケース数
	nib <- len-nia						# 第 2 群のケース数
	na <- c(nia, (rep(nia, k)-cumsum(tg[,1]+tg[,2]))[-k])	# 第 1 群の総症例数
	nb <- c(nib, (rep(nib, k)-cumsum(tg[,3]+tg[,4]))[-k])	# 第 2 群の総症例数
	m <- tg[,2]+tg[,4]					# 両群の死亡数の合計
	s <- m > 0						# 死亡があった死亡期間のベクトル
	m <- m[s]						# 取り出す
	r <- (na+nb)[s]						# 各死亡期間における死亡リスク保有者数
	A <- nb[s]/r
	U <- sum(tg[,4])-sum(m*A)				# 検定統計量
	I <- sum(m*(r-m)/(r-1)*A*(1-A))				# 誤差分散
	Z <- U/sqrt(I)						# 検定統計量を標準得点に換算
	P <- pnorm(abs(Z), lower.tail=FALSE)*2			# P 値
	return(structure(list(statistic=c(U=U, "V(U)"=I, Z=Z), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# Deming 法による回帰直線のパラメータ推定
Deming <- function(	x,						# 独立変数ベクトル
			y,						# 従属変数ベクトル
			n=1,						# ブートストラップ法で信頼区間を求めるときの回数
			a=1)						# 分散比
{
	Deming0 <- function(x, y)					# 1 組のデータについて，切片と傾きの推定値を計算する
	{
		sxx <- sum((x-mean(x))^2)
		syy <- sum((y-mean(y))^2)
		sxy <- sum((x-mean(x))*(y-mean(y)))
		if (sxy != 0) {
			Slope <- (syy-a*sxx+sqrt((syy-a*sxx)^2+4*a*sxy^2))/(2*sxy)
			Intercept <- mean(y)-Slope*mean(x)
		} else {
			Slope <- Intercept <- NA
		}
		return(c(Intercept=Intercept, Slope=Slope))
	}
	Driver <- function(x, y)					# ブートストラップ法のためのドライバー
	{
		n <- length(x)
		suffix <- sample(n, n, replace=TRUE)			# リサンプリング
		return(Deming0(x[suffix], y[suffix]))			# リサンプリングしたデータについてパラメータを推定
	}
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	OK <- complete.cases(x, y)					# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	ans <- list(coefficients=Deming0(x, y),				# 引数に対してパラメータを推定する
		    names.xy=names.xy, x=x, y=y)
	if (n > 1) {
		ans2 <- replicate(n, Driver(x, y))			# ブートストラップを n 回実行
		ans <- append(ans, list(intercepts=ans2[1,], slopes=ans2[2,]))
	}
	class(ans) <- "deming"						# print, plot メソッドのためにクラス名をつけておく
	return(ans)
}
# print メソッド
print.deming <- function(	obj,					# "deming" オブジェクト
				digits=5,				# 表示桁数
				sig=0.95)				# 信頼度
{
	if (length(obj) == 4) {
		cat("Intercept:", round(obj$coefficients[1], digits),
		    "    Slope:", round(obj$coefficients[2], digits), "\n")
	}
	else {
		alpha <- (1-sig)/2
		LCL <- c(quantile(obj$intercepts, alpha, na.rm=TRUE), quantile(obj$slopes, alpha, na.rm=TRUE))
		UCL <- c(quantile(obj$intercepts, 1-alpha, na.rm=TRUE), quantile(obj$slopes, 1-alpha, na.rm=TRUE))
		ans <- data.frame(obj$coefficients, LCL, UCL)
		dimnames(ans) <- list(c("Intercept:", "Slope:"),
				      c("Estimate", paste(c(alpha, 1-alpha), "%", sep="")))
		print(round(ans, digits=digits))
	}
}
# plot メソッド
plot.deming <- function(obj,						# "deming" オブジェクト
			posx="topleft", posy=NULL,			# legend 関数のための位置引数
			xlab=obj$names.xy[1], ylab=obj$names.xy[2],	# 軸の名前
			hist=FALSE,					# ヒストグラムを描くとき TRUE にする
			...)						# その他の任意の plot 関数の引数
{
	if (hist && length(obj) == 6) {					# ブートストラップの結果を，hist=TRUE のときに，ヒストグラムで表示する
		layout(matrix(1:2, 2))
		hist(obj$intercepts, xlab="Intercept", main="", right=FALSE)
		hist(obj$slopes, xlab="Slope", main="", right=FALSE)
		layout(1)
	}
	else {								# 散布図と Deming 法の回帰直線と直線回帰式を表示する
		plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
		abline(obj$coefficients)
		abline(lm(obj$y~obj$x), lty=2, col=2)
		legend(posx, posy, legend=c("Deming", "linear regression"), lty=1:2, col=1:2)
	}
}
# n 番目のフィボナッチ数を求める（n はベクトルでもよい）
Fibonacci <- function(n)
{
	return( (((1+sqrt(5))/2)^n - ((1-sqrt(5))/2)^n) / sqrt(5) )
}
# 対数尤度比に基づく独立性の検定を行う（htest クラスの結果を返す）
G2 <- function(	mat,						# 分割表（合計欄を除く）
		correct=FALSE)					# Williams の連続性の補正を行うときに TRUE にする
{
	ln <- function(n) sum(ifelse(n == 0, 0, n*log(n)))	# n*ln(n) を適切に行う関数
	data.name <- deparse(substitute(mat))
	method <- "対数尤度比に基づく独立性の検定（G-squared test）"
	n <- sum(mat)						# 全サンプルサイズ
	n1 <- rowSums(mat)					# 行和
	n2 <- colSums(mat)					# 列和
	G2 <- 2*(ln(mat)-ln(n1)-ln(n2)+ln(n))			# G 統計量
	a <- nrow(mat)						# 分割表の行数
	b <- ncol(mat)						# 分割表の列数
	df <- (a-1)*(b-1)					# G の自由度
	if (correct == TRUE) {					# 連続性の補正
		method <- paste(method, "連続性の補正")
		G2 <- G2/(1+(n*sum(1/n1)-1)*(n*sum(1/n2)-1)/(6*n*a*b))
	}
	P <- pchisq(G2, df, lower.tail=FALSE)			# P 値
	names(G2) <- "G-squared"
	names(df) <- "df"
	return(structure(list(statistic=G2, parameter=df, 
		p.value=P, method=method, data.name=data.name, observed=mat), 
		class="htest"))					# 結果をまとめて返す
}
# 一般化 Wilcoxon 検定を行う
Gen.Wil <- function(	group,							# 群を識別するベクトル
			event,							# 死亡なら 1，生存なら 0
			time)							# 生存期間ベクトル
{
	getU <- function(tx, sx, ty, sy)					# U の計算
	{
		if ((tx < ty && sx == 1 && sy == 1) | (tx <= ty && sx == 1 && sy == 0)) {
			return(-1)
		}
		else if ((tx > ty && sx == 1 && sy == 1) | (tx >= ty && sx == 0 && sy == 1)) {
			return(1)
		}
		else {
			return(0)
		}
	}
	method <- "一般化 Wilcoxon 検定"
	data.name <- sprintf("time: %s, event: %s, group: %s",
		deparse(substitute(time)), deparse(substitute(event)), deparse(substitute(group)))
	OK <- complete.cases(group, event, time)				# 欠損値を持つケースを除く
	group <- group[OK]
	event <- event[OK]
	time <- time[OK]
	n <- length(group)
	o <- order(group)							# グループによって並べ替え
	time <- time[o]
	group <- group[o]
	event <- event[o]
	na <- table(group)[1]							# 各群のデータ個数
	nb <- n-na								# 各群のデータ個数

	W <- 0									# 検定統計量
	for (i in 1:na) {
		for (j in (na+1):n) {
			W <- W+getU(time[i], event[i], time[j], event[j])
		}
	}

	Var.W <- 0								# 分散
	for (i in 1:n) {
		temp <- 0
		for (j in 1:n) {
			temp <- temp+getU(time[i], event[i], time[j], event[j])
		}
		Var.W <- Var.W+temp^2
	}
	Var.W <- na*nb*Var.W/n/(n-1)
	Z <- W/sqrt(Var.W)							# Z 値
	P <- pnorm(abs(Z), lower.tail=FALSE)*2					# P 値
	return(structure(list(statistic=c(W=W, "Var(W)"=Var.W, Z=Z), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# ローレンツ曲線を描き，ジニ係数を計算する
Gini.index <- function(	y,			# データベクトル
			main="",		# 図のタイトル（省略時は何も書かない）
			xlab="",		# x 軸の名前（省略時は何も書かない）
			ylab="")		# y 軸の名前（省略時は何も書かない）
{
	stopifnot(y >= 0)			# 非負データでなければならない
	n <- length(y)				# データの個数
	y <- sort(y)				# 小さい順に並べる
	y <- cumsum(y)				# 累積度数をとる
	y <- c(0, y/y[n])			# 累積相対度数（先頭に 0 を加える）
	x <- seq(0, 1, length=n+1)		# 0 ～ 1 を等間隔に区切ったベクトルを作る
	old <- par(xaxs="i", yaxs="i")
	plot(x, y, type="l", col="blue",	# これを結ぶとローレンツ曲線
		main=main, xlab=xlab, ylab=ylab)
	abline(0, 1)				# 対角線（原点を通る，傾き 1 の直線）を描く
	par(old)
	return(2*sum(x-y)/n)			# ジニ係数
}
# ホッジス・レーマン推定量（中央値の推定）
HLe <- function(x)						# データベクトル
{
	x <- x[!is.na(x)]					# 欠損値を除く
	temp <- outer(x, x, "+")/2				# データを二つずつ取り出し，平均値を求める
	return(median(temp[upper.tri(temp, diag=TRUE)]))	# その中央値を求める
}
# Shannon-Wienwer の多様度指数 H' を計算する
Heterogeneity <- function(	x,					# 観察度数ベクトル
				base=NULL)				# 対数の底（デフォルトは自然対数）
{
	x <- x[x > 0] / sum(x)						# 割合に変換
	return(-sum(x*log(x)*ifelse(is.null(base), 1, 1/log(base))))	# 結果を返す
}
# ヨンキー検定
Jonckheere <- function(	x,						# データベクトル
			g,						# 群変数ベクトル
			correct=FALSE,					# 連続性の補正
			alternative = c("increasing", "decreasing"))	# 帰無仮説の種類（必ず片側検定）
{
	method <- "ヨンキー検定"
	data.name <- paste(deparse(substitute(x)), "by", deparse(substitute(g)))
	sgn <-  ifelse(match.arg(alternative) == "increasing", "<", ">")
	alternative <- sprintf("control %s= treat-1 %s= ... %s= treat-n", sgn, sgn, sgn)
	OK <- complete.cases(x, g)					# 欠損値を持つケースを除く
	x <- x[OK]
	g <- g[OK]
	ni <- table(g)							# 各群のケース数
	gv <- as.numeric(names(ni))					# 群を表す数値
	a <- length(ni)							# 群の個数
	n <- sum(ni)							# 全ケース数
	sumSij <- sumTij <- 0
	for (i in 1:(a-1)) {
		di <- x[g==gv[i]]
		for (j in (i+1):a) {
			dj <- x[g==gv[j]]
			sumTij <- sumTij+sum(outer(di, dj, sgn))
			sumSij <- sumSij+sum(outer(di, dj, "=="))
		}
	}
	tau <- table(x)
	V <- (n*(n-1)*(2*n+5)-sum(ni*(ni-1)*(2*ni+5))-sum(tau*(tau-1)*(2*tau+5)))/72 +
		sum(ni*(ni-1)*(ni-2))*sum(tau*(tau-1)*(tau-2))/(36*n*(n-1)*(n-2)) +
		sum(ni*(ni-1))*sum(tau*(tau-1))/(8*n*(n-1))		# 分散
	J0 <- sumTij+sumSij/2						# 検定統計量
	E <- (n^2-sum(ni^2))/4						# 期待値
	Z <- (abs(J0-E)-0.5*correct)/sqrt(V)				# 標準化得点
	P <- pnorm(Z, lower.tail=FALSE)					# P 値
	return(structure(list(statistic=c(J=J0, "E(J)"=E, "V(J)"=V, "Z-value"=Z), p.value=P,
		method=method, alternative=alternative, data.name=data.name), class="htest"))
}
# Linear-by-Linear検定（Mantel の傾向検定）
Mantel <- function(	r.i,			# 反応数のベクトル
			n.i,			# ケース数のベクトル
			x.i=seq(along=r.i))	# 外的基準。省略されたときは1から始まる整数値が仮定される
{
	data.name <- paste(deparse(substitute(r.i)), "out of", deparse(substitute(n.i)),
					",\n using scores:", paste(x.i, collapse = " "))
	method <- "Linear-by-Linear検定（Mantel の傾向検定）"
	x <- rep(c(x.i, x.i), c(r.i, n.i-r.i))  # 外的基準の展開
	n.r <- sum(r.i)				# 反応総数
	y <- rep(1:2, c(n.r, sum(n.i)-n.r))	# 反応あり・なしの 2 値データを展開
	s <- (sum(n.i)-1)*cor(x, y)^2           # 検定統計量
	df <- 1
	p <- pchisq(s, df, lower.tail=FALSE)
	names(s) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=s, parameter=df, p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# Major Axis regression（主成分回帰）
MA <- function(x, y)							# 2 つの変数
{
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	OK <- complete.cases(x, y)					# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	s2 <- cov(cbind(x, y))						# 分散・共分散行列
	ev <- eigen(s2)$values						# 固有値
	b <- s2[1, 2]/(ev[1]-s2[2, 2])					# 傾き
	a <- mean(y)-b*mean(x)						# 切片
	result <- list(names.xy=names.xy, x=x, y=y, ev=ev, intercept=a, slope=b)
	class(result) <- "MA"
	return(result)
}
# print メソッド
print.MA <- function(	obj,						# "MA" オブジェクト
			digits=5,					# 表示桁数
			sig=0.95)					# 信頼度
{
	n <- length(obj$x)						# ケース数
	ev <- obj$ev
	H <- qf(sig, 1, n-2)/(ev[1]/ev[2]+ev[2]/ev[1]-2)/(n-2)
	CL <- tan(atan(obj$slope)+c(-0.5, 0.5)*asin(2*sqrt(H)))			# 傾きの信頼限界値
	ans <- data.frame(c(obj$intercept, obj$slope),
			  c(NA, CL[1]), c(NA, CL[2]))
	alpha <- 0.5-sig/2
	dimnames(ans) <- list(c("Intercept:", "Slope:"),
			      c("Estimate", paste(c(alpha, 1-alpha), "%", sep="")))
	print(ans, digits=digits)		
}
# plot メソッド
plot.MA <- function(	obj,						# "MA" オブジェクト
			posx="topleft", posy=NULL,			# legend 関数のための位置引数
			xlab=obj$names.xy[1],				# x 軸の名前
			ylab=obj$names.xy[2],				# y 軸の名前
			...)						# その他の任意の plot 関数の引数

{
	plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
	abline(obj$intercept, obj$slope)
	abline(lm(obj$y~obj$x), lty=2, col=2)
	legend(posx, posy, legend=c("Major Axis", "linear regression"), lty=1:2, col=1:2)
}
# マクネマー検定
McNemar <- function(tbl)			# 分割表
{
	data.name <- deparse(substitute(tbl))
	method <- "拡張されたマクネマー検定（二項検定に帰着）"
	n1 <- sum(tbl[lower.tri(tbl)])		# x の方が大きいデータ対の数
	n2 <- sum(tbl[upper.tri(tbl)])		# y の方が大きいデータ対の数
	p <- binom.test(c(n1, n2))$p.value	# 二項検定の両側確率
	names(n1) <- "n1"
	names(n2) <- "n2"
	return(structure(list(statistic=n1, parameter=n2, p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# Morisita の類似度指数 Cλ を計算する
Morisita <- function(	a,			# 第一群の観察度数ベクトル
			b)			# 第二群の観察度数ベクトル
{
	stopifnot(length(a) == length(b))
	na <- sum(a)
	nb <- sum(b)
	return(2*sum(a*b)/((sum(a*(a-1)/na/(na-1))+sum(b*(b-1)/nb/(nb-1)))*na*nb))
}
# 負の超幾何分布

# ガンマ関数の自然対数を使って階乗を計算
# n! = exp(lgamma(n+1))
NegativeGeometric <- function(x, N, n, r)
{
	exp(lgamma(x)-lgamma(r)-lgamma(x-r+1)-lgamma(N+1)+lgamma(n+1)+lgamma(N-n+1)+lgamma(N-x+1)-lgamma(n-r+1)-lgamma(N-x-n+r+1))
}

# 二項係数 choose(n, k) の自然対数を利用
# nCk = exp(lchoose(n, k))

NegativeGeometric2 <- function(x, N, n, r)
{
	exp(lchoose(x-1, r-1)+lchoose(N-x, n-r)-lchoose(N, n))
}

# パレート図を描く
Pareto <- function(	x,		# 度数分布表
			ymax=sum(x),	# 度数軸の最大値
			sort.flag=TRUE,	# 度数分布の大きい順に並べ替えるかどうか
			col=NULL,	# 度数を表す矩形の描画色
			density=NULL,	# 度数を表す矩形のハッチ
			lwd=1,		# 累積度数曲線の描画線種
			las=0,		# 軸のラベルの描き方
			main=NULL,	# グラフのタイトル
			xlab=NULL,	# 横軸のラベル
			ylab="度数",	# 度数軸のラベル
			ylab2="累積％")	# 累積度数軸のラベル
{
	if (sort.flag) {
		x <- sort(x, decreasing=TRUE)
	}
	old <- par(xpd=TRUE, mar=c(5, 5, 1, 5)+0.1)
	barplot(x, space=0, xlab=xlab, ylab=ylab, ylim=c(0, ymax),
	        las=las, col=col, density=density, main=main)
	lines(cumsum(x), lwd=lwd)
	par(old)
	axis(4, at=seq(0, sum(x), length=11), labels=0:10*10,
	     pos=length(x)+0.4, las=las)
	mtext(ylab2, side=4, las=las)
}
# 線形判別関数における説明変数の有意性検定
PartialF <- function(	data,					# 説明変数のデータ行列（データフレーム）
			group,					# 各行がどの群に属するかを表す変数
			use)					# 分析に使用する変数番号のベクトル
{
	data <- data[, use, drop=FALSE]
	p <- ncol(data)
	ng <- nlevels(group)
	ncase <- length(group)
	split.data <- split(data, group)
	num <- sapply(split.data, nrow)
	t <- var(data)*(ncase-1)
	w <- matrix(colSums(t(matrix(sapply(split.data, var), ncol=ng))*(num-1)), p)
	f <- diag(solve(t))/diag(solve(w))
	ni <- length(use)
	idf1 <- ng-1
	idf2 <- ncase-(ng-1)-ni
	f <- idf2 / idf1 * (1-f) / f
	p <- pf(f, idf1, idf2, lower.tail=FALSE)
	result <- cbind(f, p)
	dimnames(result) <- list(colnames(data), c("Partial F", "p-value"))
	return(result)
}
# Passing & Bablok 法による回帰直線のパラメータ推定
PassingBablok <- function(	x, # 独立変数ベクトル
						y, # 従属変数ベクトル
						n=1) # ブートストラップ法で信頼区間を求めるときの回数
{
  PassingBablok0 <- function(x, y) # 1 組のデータについて，切片と傾きの推定値を計算する
  {
    cor.sign2 = sign(cor(x, y, method="kendall"))
    suffix <- combn(length(x), 2)
    slope <- diff(matrix(y[suffix], 2))/diff(matrix(x[suffix], 2)) # すべての二点の組合せの傾きを計算
    slope <- slope[!is.nan(slope)] # 傾きが垂直（計算上は NaN になる）のものを除く
    if (cor.sign2 == cor.sign) { # Passing & Bablok 法では傾きが -1 になるものを除き，
      slope <- slope[slope != -1]
      k <- sum(slope < -1) # 傾きが -1 未満のものの個数が必要
      slope <- sort(slope) # メディアンを推定値にするが，ちょっと特殊なことをする
    }
    else {
      slope <- slope[slope != 1]
      k <- sum(slope > 1)
      slope <- sort(slope, decreasing=TRUE) # メディアンを推定値にするが，ちょっと特殊なことをする
    }
    n <- length(slope)
    if (n %% 2 == 0) { # 普通のメディアンの計算より，k 番目に大きいものを使う
      Slope <- (slope[(n+1)/2+k]+slope[(n+1)/2+k+1])/2
    }
    else {
      Slope <- slope[(n+1)/2+k]
    }
    return(c(Intercept=median(y-Slope*x), Slope=Slope)) # 結果を返す（2 要素を持つベクトル）
  }
  Driver <- function(x, y) # ブートストラップ法のためのドライバー
  {
    n <- length(x)
    suffix <- sample(n, n, replace=TRUE) # リサンプリング
    return(PassingBablok0(x[suffix], y[suffix])) # リサンプリングしたデータについてパラメータを推定
  }
  names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
  OK <- complete.cases(x, y) # 欠損値を持つケースを除く
  x <- x[OK]
  y <- y[OK]
  cor.sign = sign(cor(x, y, method="kendall"))
  ans <- list(coefficients=PassingBablok0(x, y), # 引数に対してパラメータを推定する
			names.xy=names.xy, x=x, y=y)
  if (n > 1) {
    ans2 <- replicate(n, Driver(x, y)) # ブートストラップを n 回実行
    ans <- append(ans, list(intercepts=ans2[1,],  slopes=ans2[2,]))
  }
  class(ans) <- "passing.bablok" # print, plot メソッドのためにクラス名をつけておく
  return(ans)
}
# print メソッド
print.passing.bablok <- function(	obj, # "passing.bablok" オブジェクト
								digits=5, # 表示桁数
								sig=0.95) # 信頼度
{
  if (length(obj) == 4) {
    cat("Intercept:", round(obj$coefficients[1], digits),
	   "    Slope:", round(obj$coefficients[2], digits), "\n")
  }
  else {
    alpha <- (1-sig)/2
    LCL <- c(quantile(obj$intercepts, alpha, na.rm=TRUE), quantile(obj$slopes, alpha, na.rm=TRUE))
    UCL <- c(quantile(obj$intercepts, 1-alpha, na.rm=TRUE), quantile(obj$slopes, 1-alpha, na.rm=TRUE))
    ans <- data.frame(obj$coefficients, LCL, UCL)
    dimnames(ans) <- list(c("Intercept:", "Slope:"),
					   c("Estimate", paste(c(alpha, 1-alpha), "%", sep="")))
    print(ans, digits=digits)		
  }
}
# plot メソッド
plot.passing.bablok <- function(obj, # "passing.bablok" オブジェクト
							posx="topleft", posy=NULL, # legend 関数のための位置引数
							xlab=obj$names.xy[1], # x 軸の名前
							ylab=obj$names.xy[2], # y 軸の名前
							hist=FALSE, # ヒストグラムを描くとき TRUE にする
							...) # その他の任意の plot 関数の引数
{
  if (hist && length(obj) == 6) { # ブートストラップの結果を，hist=TRUE のときに，ヒストグラムで表示する
    layout(matrix(1:2, 2))
    hist(obj$intercepts, xlab="Intercept", main="", right=FALSE)
    hist(obj$slopes, xlab="Slope", main="", right=FALSE)
    layout(1)
  }
  else { # 散布図と Passing & Bablok 法の回帰直線と直線回帰式を表示する
    plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
    abline(obj$coefficients)
    abline(lm(obj$y~obj$x), lty=2, col=2)
    legend(posx, posy, legend=c("Passing-Bablok", "linear regression"), lty=1:2, col=1:2)
  }
}
# ペリの方法
Peritz <- function(	data,						# データベクトル
			group,						# 群別ベクトル
			alpha=0.05,					# 有意水準
			statistics=c("F", "Q"))				# 検定統計量の種類
{
# 書式付き print 関数
	printf <- function(fmt, ...)					# 書式と項目
	{
		cat(sprintf(fmt, ...))
	}

# 一つずつの検定結果
	judge <- function(	P,					# 検定対象となる群番号を表す数字の集合
				p,					# 検定対象となる群の個数（P の要素数）
				stat.P)					# 検定統計量（Q.P または F.P）
	{
		judge2 <- function(P, ns.list, stat.P, c.p)		# 下請け関数
		{
			if (!is.null(ns.list) &&			# 保留された仮説のリストがあって，
		   	 any(mapply(function(tbl)			# P がそれらの仮説の部分集合であるときは，
					all(is.element(P, tbl)), ns.list))) {
				return("NSI")				# この仮説を誘導する仮説が保留されているのでこの帰無仮説も保留する
			}
			else if (stat.P >= c.p) {			# 検定統計量が棄却限界値より大きければ，
				return("S")				# 帰無仮説を棄却する
			}
			else {						# 棄却限界値未満ならば，
				return("NS")				# 帰無仮説を保留する
			}
		}

		res <- judge2(P, NS.list, stat.P, c.p[p])		# まずは，テューキー・ウェルシュの方法により判定
		if (res == "NS") {					# 帰無仮説を保留するときは，
			NS.list <<- append(NS.list, list(P))		# 群番号の集合を保留された仮説のリストに付け加える
		}
		res.NK <- judge2(P, NS.list.NK, stat.P, c.p.NK[p])	# 次に，ニューマン・コイルスの方法により判定
		if (res.NK == "NS") {					# 帰無仮説を保留するときは，
			NS.list.NK <<- append(NS.list.NK, list(P))	# 群番号の集合を保留された仮説のリストに付け加える
		}
		n.res <<- n.res+1					# 結果を蓄積する
		result[n.res,] <<- c(paste(P, collapse=","), p, stat.P, alpha.p[p], c.p[p], res, alpha, c.p.NK[p], res.NK)
	}

# Q 統計量を使う検定関数
	funcQ <- function(P)						# 群番号を表す数字の集合
	{
		p <- length(P)						# 群番号を表す数字の集合 P の要素数
		pair <- combn(p, 2, function(i) P[i])			# P から 2 つの群を取り出す組み合わせ
		Q.P <- max(apply(pair, 2,
				 function(ij) QP[ij[2], ij[1]]))	# 検定統計量（事前に Q.P として計算してある）
		judge(P, p, Q.P)					# 検定結果の評価
	}

# F 統計量を使う検定関数
	funcF <- function(P)						# 群番号を表す数字の集合
	{
		p <- length(P)						# 群番号を表す数字の集合 P の要素数
		np <- n[P]						# P で表される群のサンプルサイズ
		mp <- m[P]						# P で表される群の平均値
		nt <- sum(np)						# 全体のサンプルサイズ
		mean.p <- sum(np*mp)/nt					# 全体の平均値
		vb.p <- sum(np*(mp-mean.p)^2)				# 群間分散
		F.P <- vb.p/(p-1)/vw					# F 統計量
		judge(P, p, F.P)					# 検定結果の評価
	}

# テューキー・ウェルシュとニューマン・コイルスの組み合わせ
	proc11 <- function(i)						# 食い違う結果があるデータフレームの行番号
	{
		P <- as.integer(unlist(strsplit(result[i,1], split=","))) # P
		p <- length(P)
		if (!is.null(NS.list.Peritz) &&				# 保留された仮説のリストがあって，
		    any(mapply(function(tbl)				# P がそれらの仮説の部分集合であるときは，
				all(is.element(P, tbl)), NS.list.Peritz))) {
			return("NS")					# この仮説を誘導する仮説が保留されているのでこの帰無仮説も保留する
		}
		else if (result[i, 3] >= result[i, 5]) {		# ZP ≧ cp なら帰無仮説を棄却する
			return("S")
		}
		else {
			comprimentP <- (1:a)[-P]			# P の補集合
			m <- length(comprimentP)			# P の補集合の要素数
			sufix <- NULL					# チェックすべき結果のある行番号
			for (i in m:2) {
				temp <- combn(m, i,			# P の補集合の全ての部分集合について，
					      function(arg) paste(comprimentP[arg], collapse=","))
				sufix <- c(sufix, sapply(temp,		# 文字列に一致する行番号を探す
					      function(temp) which(result[,1]==temp)))
			}
			temp <- result[sufix, 3] >= result[sufix, 5]	# ZP ≧ cp か？
			if (all(temp)) {				# 全てが真なら，
				return("S")				# 帰無仮説を棄却する
			}
		}
		NS.list.Peritz <<- append(NS.list.Peritz, list(P))	# 群番号の集合を保留された仮説のリストに付け加える
		return("NS")
	}

# 結果を出力する
	print.out <- function(result)					# データフレームに納めた結果
	{
		for (i in 1:nrow(result)) {
			strP <- paste("{", result[i,1], "}", sep="")
			if (i==1 || result[i, 2] != result[i-1, 2]) {	# p が同じである P の最初の検定結果
				printf(format1, strP, result[i,2], result[i,3], result[i,4], result[i,5],
					result[i,6], result[i,7], result[i,8], result[i,9], result[i,10])
			}
			else {						# 二番目以降の検定結果
				printf(format2, strP, result[i,3], result[i,6], result[i,9], result[i,10])
			}
		}
		result2 <- result[result[,2] == 2,c(1,10)]		# 二群比較の部分
		result3 <- matrix("--- ", a, a)				# 結果表
		sapply(1:nrow(result2),
				function(i) {
					temp <- ifelse(result2[i,2] == "S", "*   ", "n.s.")
					ij <- as.integer(unlist(strsplit(result2[i,1], split=",")))
					result3[ij[1], ij[2]] <<- result3[ij[2], ij[1]] <<- temp 
				})
		colnames(result3) <- rownames(result3) <- paste("A", 1:a, sep="")
		printf("\n判定結果\n")
		print(result3, quote=FALSE)
	}

# 関数本体
	statistics <- match.arg(statistics)				# 引数の補完
	OK <- complete.cases(data, group)				# 欠損値を持つケースを除く
	data <- data[OK]
	group <- group[OK]
	n <- tapply(data, group, length)
	m <- tapply(data, group, mean)
	s <- tapply(data, group, sd)
	a <- length(n)							# 群の数
	nt <- sum(n)							# 全体の標本サイズ
	dfw <- nt-a							# 群内平方和の自由度
	vw <- sum((n-1)*s^2)/dfw					# 群内分散
	NS.list <- NULL							# 検定結果を保留した仮説のリスト
	NS.list.NK <- NULL						# 検定結果を保留した仮説のリスト
	NS.list.Peritz <- NULL						# 検定結果を保留した仮説のリスト
	result <- matrix("", 2^a-a-1, 9)
	n.res <- 0
	format0 <- sprintf("%%-%is %%s %%6s %%7s %%7s %%-4s %%7s %%7s %%-4s  %%-6s\n", 2*a+1)		# 結果の出力書式0
	format1 <- sprintf("%%-%is %%i %%6.3f %%7.4f %%7.3f %%-4s %%7.4f %%7.3f %%-4s  %%-6s\n", 2*a+1)	# 結果の出力書式1
	format2 <- sprintf("%%-%is   %%6.3f                 %%-4s                 %%-4s  %%-6s\n", 2*a+1)	# 結果の出力書式2
	alpha.p <- c(NA, 1-(1-alpha)^(2:(a-2)/a), alpha, alpha)		# αp を計算
	if (a == 3) alpha.p <- c(NA, alpha, alpha)
	if (statistics == "Q") {
		QP <- matrix(NA, a, a)					# 二群の検定統計量を前もって計算しておく領域を確保
		QP[lower.tri(QP)] <- combn(a, 2, function(ij) {
			i <- ij[1]
			j <- ij[2]
			abs(m[i]-m[j])/sqrt(vw*(1/n[i]+1/n[j]))		# 検定統計量
		})
		qTukey <- qtukey(alpha.p[1:a], 1:a, dfw, lower.tail=FALSE)
		c.p <- c(NA, qTukey[2]/sqrt(2))				# cp を計算
		sapply(3:a, function(i) c.p[i] <<- max(qTukey[i]/sqrt(2), c.p[i-1]))
		qTukey.NK <- c(NA, qtukey(alpha, 2:a, dfw, lower.tail=FALSE))
		c.p.NK <- c(NA, qTukey.NK[2]/sqrt(2))			# cp.NK を計算
		sapply(3:a, function(i) c.p.NK[i] <<- max(qTukey.NK[i]/sqrt(2), c.p.NK[i-1]))
		printf(format0, "P", "p", "QP", "αp", "cp(TW)", "判定", "αp", "cp(NK)", "判定", "ペリの判定")
		sapply(a:2, function(i) combn(a, i, funcQ))		# a 群から a～2 群を取り出す組み合わせについて検定を実行
	}
	else { # if (statistics == "F") {
		c.p <- qf(alpha.p, 0:(a-1), dfw, lower.tail=FALSE)
		c.p.NK <- c(NA, qf(alpha, 1:(a-1), dfw, lower.tail=FALSE))
		printf(format0, "P", "p", "FP", "αp", "cp(TW)", "判定", "αp", "cp(NK)", "判定", "ペリの判定")
		sapply(a:2, function(i) combn(a, i, funcF))		# a 群から a～2 群を取り出す組み合わせについて検定を実行
	}
	result <- data.frame(result)
	old <- options(width=256)
	con <- textConnection(capture.output(result))
	result <- read.table(con)
	close(con)
	options(old)
	result[,1] <- as.character(result[,1])
	result[,6] <- as.character(result[,6])
	result[,9] <- as.character(result[,9])
	result[,10] <- ifelse((result[,6] == "NS" | result[,6] == "NSI") & result[,9] == "S", "不一致", ifelse(result[,6] == "NSI", "NS", result[,6]))
	colnames(result) <- c("P", "p", "ZP", "αp", "cp.TW", "判定.TW", "αp.NK", "cp.NK", "判定.NK", "判定.Peritz")
	sapply(1:nrow(result), function(i) {
		if (result[i, 10] == "不一致") {
			result[i, 10] <<- proc11(i)
		}
	})
	print.out(result)
	invisible(result)
}
# ポリヤ・エッゲンベルガー分布の確率を計算する
# n, p, δ を与える場合
PolyaEggenberger <- function(	x,			# 確率変数
				n,			# 標本サイズ
				p,			# 母比率
				delta)			# 追加する割合 δ
{
	exp(						# 対数で計算して最後に逆対数をとる
		sum(lchoose(n, x),
		sapply(0:x, function(i) ifelse(i == 0, 0, log(p+(i-1)*delta)-log(1+(i-1)*delta))),
		sapply((x+1):n, function(i) log(1-p+(i-x-1)*delta)-log(1+(i-1)*delta)))
	)
}
# λ, r を与える場合
PolyaEggenberger2 <- function(	x,			# 確率変数
				lambda,			# λ = n*p
				r)			# r = n*δ
{
	exp(						# 対数で計算して最後に逆対数をとる
		sum(sapply(0:(x-1), function(i) ifelse(i < 0, 0, log(1+i*r/lambda))))
		+x*log(lambda)
		-lgamma(x+1)
		+(-lambda/r-x)*log(1+r)
	)
}
# ポリヤ・エッゲンベルガー分布への適合度の検定
PolyaEggenbergerdist <- function(d,             # 度数ベクトル
                                 x=NULL)        # 階級値ベクトル
{
  PolyaEggenberger2 <- function(x,              # 確率変数
                                lambda,         # λ = n*p
                                r)              # r = n*δ
  {
    exp(                                        # 対数で計算して最後に逆対数をとる
      sum(sapply(0:(x-1), function(i) ifelse(i < 0, 0, log(1+i*r/lambda))))
      +x*log(lambda)
      -lgamma(x+1)
      +(-lambda/r-x)*log(1+r)
    )
  }
#===
  if (is.null(x)) {
    stop("関数の仕様が変更されました。度数ベクトルと同じ長さの階級値ベクトルも指定してください。")
  }
  if (length(x) != length(d)) {
    stop("度数ベクトル階級値ベクトルの長さは同じでなければなりません。")
  }
  data.name <- paste(deparse(substitute(d)), deparse(substitute(x)), sep=", ")
  method <- "ポリヤ・エッゲンベルガー分布への適合度の検定"

  o <- numeric(diff(range(x))+1)
  o[x-min(x)+1] <- d
  x <- min(x):max(x)

  k <- length(o)                                # 階級数
  n <- sum(o)                                   # データ数
  lambda <- sum(o*x)/n                          # 平均値（=λ）
  r <- sum(o*(x-lambda)^2)/n/lambda-1           # パラメータ r
  p <- sapply(x, PolyaEggenberger2, lambda, r)  # 確率
  if (min(x) != 0) {                            # 最初と最後の階級値の確率は階級値以下・以上の確率を併合する
    p[1] = sum(sapply(0:min(x), PolyaEggenberger2, lambda, r))
  }
  p[k] <- 1-sum(p[-k])                          # 最後の確率は合計が 1 になるように
  e <- n*p                                      # 期待値
  table <- data.frame(x, o, p, e)               # 結果をデータフレームにまとめる
  rownames(table) <- paste("c-", x, sep="")

  while (e[1] < 1) {                            # 1 未満のカテゴリーの併合
    o[2] <- o[2]+o[1]
    e[2] <- e[2]+e[1]
    o <- o[-1]
    e <- e[-1]
    k <- k-1
  }
  while (e[k] < 1) {                            # 1 未満のカテゴリーの併合
    o[k-1] <- o[k-1]+o[k]
    e[k-1] <- e[k-1]+e[k]
    o <- o[-k]
    e <- e[-k]
    k <- k-1
  }
  chisq <- sum((o-e)^2/e)                       # カイ二乗統計量
  df <- k-3                                     # 自由度
  p <- pchisq(chisq, df, lower.tail=FALSE)      # P 値
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, parameter=df, p.value=p,
    estimate=c(n=n, lambda=lambda, r=r), method=method,
    data.name=data.name, table=table), class=c("htest", "PolyaEggenbergerdist")))
}
# summary メソッド
summary.PolyaEggenbergerdist <- function(obj,   # PolyaEggenbergerdist が返すオブジェクト
                                         digits=5)
{
  table <- obj$table
  colnames(table) <- c("階級", "度数", "確率", "期待値")
  cat("\n適合度\n\n")
  print(table, digits=digits, row.names=FALSE)
}
# plot メソッド
plot.PolyaEggenbergerdist <- function(obj,      # PolyaEggenbergerdist が返すオブジェクト
                                      ...)      # barplot へ渡す引数
{
  table <- obj$table
  nr <- nrow(table)
  pos <- barplot(table$o, space=0, ...)         # 観察度数を barplot で描く
  old <- par(xpd=TRUE)
  points(pos, table$e, pch=3)                   # 理論度数を，記号 + で示す
  text(pos, -strheight("H"), table$x)
  par(old)
}
# Reduced Major Axis regression
RMA <- function(x, y)							# 2 つの変数
{
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	OK <- complete.cases(x, y)					# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	n <- length(x)							# ケース数
	n1 <- n-1
	df <- n-2							# 信頼限界を求めるときに必要な t 分布の自由度
	slope <- sign(cor(x, y))*sd(y)/sd(x)				# 傾き
	intercept <- mean(y)-slope*mean(x)				# 切片
	MSE <- (var(y)-cov(x, y)^2/var(x))*n1/df			# 標準誤差
	SE.intercept <- sqrt(MSE*(1/n+mean(x)^2/var(x)/n1))		# 切片の標準誤差
	SE.slope <- sqrt(MSE/var(x)/n1)					# 傾きの標準誤差
	result <- list(names.xy=names.xy, x=x, y=y,
		       intercept=intercept, SE.intercept=SE.intercept,
		       slope=slope, SE.slope=SE.slope)
	class(result) <- "RMA"
	return(result)
}
# print メソッド
print.RMA <- function(	obj,						# "RMA" オブジェクト
			digits=5,					# 表示桁数
			sig=0.95)					# 信頼度
{
	alpha <- (1-sig)/2
	df <- length(obj$x)-2
	intercept <- obj$intercep
	slope <- obj$slope
	SE.intercept <- obj$SE.intercept
	SE.slope <- obj$SE.slope
	CLintercept <- intercept+qt(c(alpha, 1-alpha), df)*SE.intercept	# 切片の信頼限界値
	CLslope <- slope+qt(c(alpha, 1-alpha), df)*SE.slope		# 傾きの信頼限界値
	ans <- data.frame(c(intercept, slope), c(SE.intercept, SE.slope),
			  c(CLintercept[1], CLslope[1]), c(CLintercept[2], CLslope[2]))
	dimnames(ans) <- list(c("Intercept:", "Slope:"),
			      c("Estimate", "S.E.", paste(c(alpha, 1-alpha), "%", sep="")))
	print(ans, digits=digits)		
}
# plot メソッド
plot.RMA <- function(	obj,						# "RMA" オブジェクト
			posx="topleft", posy=NULL,			# legend 関数のための位置引数
			xlab=obj$names.xy[1],				# x 軸の名前
			ylab=obj$names.xy[2],				# y 軸の名前
			...)						# その他の任意の plot 関数の引数

{
	plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
	abline(obj$intercept, obj$slope)
	abline(lm(obj$y~obj$x), lty=2, col=2)
	legend(posx, posy, legend=c("Reduced Major Axis", "linear regression"), lty=1:2, col=1:2)
}
# 生データまたは度数分布表データに基づいて，ROC 曲線を描く。また，ROC 曲線下面積を計算する
# 生データがあるとき（後述する ROC 関数も必要なので注意）
ROC0 <- function(	disease,				# 疾病群の測定値ベクトル
			normal,					# 健康者群の測定値ベクトル
			lowest=NULL,				# データの最小値より小さいキリのよい数値
			width=NULL)				# 度数分布を作成する階級幅のキリのよい数値
{
	my.hist <- function(x, brks)				# R の hist 関数は，級限界の扱いがイヤラシイので自前の関数を書く
	{
		k <- length(brks)
		freq <- numeric(k)
		for (i in 1:(k-1)) {
			freq[i] <- sum(brks[i] <= x & x < brks[i+1])
		}
		freq[k] <- sum(x >= brks[k])
		freq
	}

	x <- c(disease, normal)					# データをプールする
	min.x <- min(x)						# 最小値
	max.x <- max(x)						# 最小値
	cat("最小値 x = ", min.x, "\n")
	cat("最大値 x = ", max.x, "\n\n")
	if (is.null(lowest) || is.null(width)) {
		temp <- pretty(c(disease, normal), n=min(length(disease)+length(normal), 50))
		lowest <- temp[1]
		width <- diff(temp)[1]
		cat("最小値より小さいキリのよい数値 = ", lowest, "\n")
		cat("度数分布を作成する階級幅の切りのよい数値 = ", width, "\n\n")
	}
	
	brks <- seq(lowest, max.x+width, by=width)
	ROC(brks, my.hist(disease, brks), my.hist(normal, brks))
	
}
# 度数分布表しかないとき（生データから計算されるときにも，下請けとして使う）	
ROC <- function(	x,					# 分割表の下限値のベクトル
			disease,				# 疾病群の度数分布ベクトル
			normal)					# 健康者群の度数分布ベクトル
{
	k <- length(x)						# 度数分布表の長さ
	stopifnot(k == length(disease) && k == length(normal))	# データのチェック
	Sensitivity <- c(rev(cumsum(rev(disease)))/sum(disease), 0)
	False.Positive.Rate <- c(rev(cumsum(rev(normal)))/sum(normal), 0)
	Specificity <- 1-False.Positive.Rate
	plot(False.Positive.Rate, Sensitivity, type="b")
	abline(h=c(0, 1), v=c(0, 1))
	c.index <- sum(sapply(1:k, function(i) (False.Positive.Rate[i]-False.Positive.Rate[i+1])*(Sensitivity[i+1]+Sensitivity[i])/2)) # area under ROC curve
	result <- cbind(x, disease, normal, Sensitivity[-k-1], Specificity[-k-1], False.Positive.Rate[-k-1])
	rownames(result) <- as.character(1:k)
	colnames(result) <- c("Value", "Disease", "Normal", "Sensitivity", "Specificity", "F.P. rate")
	list(result=result, c.index=c.index)
}
# リジット分析を行う
Ridit <- function(	na,					# 第一群の度数
			nb,					# 第二群の度数
			combine=FALSE)				# 両群をプールして基準にするとき TRUE
{
	ridit <- function(na)					# リジットを計算する関数
	{
		cna <- cumsum(na)
		return((na/2+c(0, cna[-length(cna)]))/sum(na))
	}

	data.name <- paste(deparse(substitute(na)), "and", deparse(substitute(nb)))
	if (combine == FALSE) {					# 片方の群を基準にする場合
		method <- "リジット分析（片方の群を基準にする場合）"
		n2 <- sum(nb)					# 比較される群 nb の総数
		r <- sum(nb*ridit(na))/n2
		Z <- abs(r-0.5)*sqrt(12*n2)
		P <- pnorm(Z, lower.tail=FALSE)*2
		names(Z) <- "Z-value"
		names(r) <- "ridit"
		return(structure(list(statistic=Z, p.value=P, estimate=r,
			method=method, data.name=data.name), class="htest"))
	}
	else {							# 群をプールして基準とする場合
		method <- "リジット分析（群をプールして基準とする場合）"
		r <- ridit(na+nb)
		n1 <- sum(na)
		ra <- sum(r*na)/n1
		n2 <- sum(nb)
		rb <- sum(r*nb)/n2
		Z <- abs(ra-rb)/sqrt((n1+n2)/(12*n1*n2))
		P <- pnorm(Z, lower.tail=FALSE)*2
		names(Z) <- "Z-value"
		return(structure(list(statistic=Z, p.value=P,
			estimate=c("ridit-a"=ra, "ridit-b"=rb),
			method=method, data.name=data.name), class="htest"))
	}
}
# 二要因の分散分析（SAB タイプ；RBFpq デザイン；被検者内計画）を行う
SAB <- function(data)						# 3次元配列のデータ
{
	N <- dim(data)[3]					# サンプルサイズ
	Na <- dim(data)[2]					# 要因 A の水準数
	Nb <- dim(data)[1]					# 要因 B の水準数
	Xbar <- as.numeric(apply(data, 1:2, mean))		# 水準の組み合わせごとの平均値
	SD <- as.numeric(apply(data, 1:2, sd)) *sqrt((N-1)/N)	# 水準の組み合わせごとの標準偏差
	v1 <- mean(data)					# 全体の平均値
	v2 <- sum((apply(data, 2, mean)-v1)^2)*Nb*N		# 要因 A の水準による変動(A)
	v3 <- sum((apply(data, 1, mean)-v1)^2)*Na*N		# 要因 B の水準による変動(B)
	v4 <- sum((apply(data, 1:2, mean)-v1)^2)*N		# 要因の水準の組み合わせによる変動
	v4.2 <- v4-v2-v3					# 要因 A, B の交互作用による変動(AxB)
	v5 <- sum(SD^2)*N					# 偶然誤差
	v6 <- sum(data)^2/(N*Na*Nb)				# 修正項
	v6.2 <- sum(apply(data, 3, sum)^2)/(Na*Nb)-v6		# 個人差による変動(S)
	v7 <- sum(apply(data, 2:3, sum)^2)/Nb-v6-v6.2-v2	# 要因 A の誤差変動(SxA)
	v8 <- sum(apply(data, c(1, 3), sum)^2)/Na-v6-v6.2-v3	# 要因 B の誤差変動(SxB)
	v9 <- v5-v6.2-v7-v8					# 交互作用による変動
	SS <- c(v6.2, v2, v7, v3, v8, v4.2, v9)
	dfs <- N-1
	dfa <- Na-1
	dfb <- Nb-1
	df <- c(dfs, dfa, dfs*dfa, dfb, dfs*dfb, dfa*dfb, dfs*dfa*dfb)
	MS <- SS/df
	P <- F <- rep(NA, 7)
	F[c(2, 4, 6)] <- MS[c(2, 4, 6)]/MS[c(3, 5, 7)]
	P[c(2, 4, 6)] <- pf(F[c(2, 4, 6)], df[c(2, 4, 6)], df[c(3, 5, 7)], lower.tail=FALSE)
	result <- data.frame(SS=SS, df=df, MS=MS, F=F, P=P)
	colnames(result) <- c("SS", "d.f.", "MS", "F value", "P value")
	rownames(result) <- c("S", "A", "SxA", "B", "SxB", "AxB", "SxAxB")
	class(result) <- c("anova.table", "data.frame")
	return(result)
}
SABC <- function(tbl) {
    d <- dim(tbl)
    n <- d[4]
    r <- d[1]
    q <- d[2]
    p <- d[3]
    X <- sum(tbl)^2/n/p/q/r
    A <- sum(apply(tbl, 3, sum)^2)/n/q/r
    B <- sum(apply(tbl, 2, sum)^2)/n/p/r
    C <- sum(apply(tbl, 1, sum)^2)/n/p/q
    S <- sum(apply(tbl, 4, sum)^2)/p/q/r
    AB <- sum(apply(tbl, c(3, 2), sum)^2)/n/r
    AC <- sum(apply(tbl, c(3, 1), sum)^2)/n/q
    BC <- sum(apply(tbl, c(2, 1), sum)^2)/n/p
    AS <- sum(apply(tbl, c(4, 3), sum)^2)/q/r
    BS <- sum(apply(tbl, c(4, 2), sum)^2)/p/r
    CS <- sum(apply(tbl, c(4, 1), sum)^2)/p/q
    ABS <- sum(apply(tbl, c(4, 2, 3), sum)^2)/r
    ACS <- sum(apply(tbl, c(4, 1, 3), sum)^2)/q
    BCS <- sum(apply(tbl, c(4, 2, 1), sum)^2)/p
    ABC <- sum(apply(tbl, c(3, 1, 2), sum)^2)/n
    ABCS <- sum(tbl^2)
    SS.A <- A - X
    SS.B <- B - X
    SS.C <- C - X
    SS.AB <- AB - A - B + X
    SS.AC <- AC - A - C + X
    SS.BC <- BC - B - C + X
    SS.ABC <- ABC - AB - AC - BC + A + B + C - X
    SS.T <- ABCS - X
    SS.S <- S - X
    SS.AS <- AS - A - S + X
    SS.BS <- BS - B - S + X
    SS.CS <- CS - C - S + X
    SS.ABS <- ABS - AB - AS - BS + A + B + S - X
    SS.ACS <- ACS - AC - AS - CS + A + C + S - X
    SS.BCS <- BCS - BC - BS - CS + B + C + S - X
    SS.ABCS <- ABCS - ABC - ABS - ACS - BCS +
      AB + AC + BC + AS + BS + CS - A - B - C - S + X
    SS <- c(SS.S, SS.A, SS.AS, SS.B, SS.BS, SS.C, SS.CS, SS.AB,
      SS.ABS, SS.AC, SS.ACS, SS.BC, SS.BCS, SS.ABC, SS.ABCS, SS.T)
    n1 <- n - 1
    p1 <- p - 1
    q1 <- q - 1
    r1 <- r - 1
    df <- c(n1, p1, p1 * n1, q1, q1 * n1, r1, r1 * n1, p1 * q1,
      p1 * q1 * n1, p1 * r1, p1 * r1 * n1, q1 * r1, q1 * r1 * n1,
      p1 * q1 * r1, p1 * q1 * r1 * n1, n * p * q * r - 1)
    MS <- SS/df
    suf <- 1:7 * 2
    p.value <- F <- rep(NA, 16)
    F.value[suf] <- MS[suf]/MS[suf + 1]
    p.value[suf] <- pf(F.value[suf], df[suf], df[suf + 1], lower.tail = FALSE)
    ANOVA.table <- data.frame(SS, df, MS, F.value, p.value)
    rownames(ANOVA.table) <- c("S", "A", "AxS", "B", "BxS", "C", "CxS",
      "AxB", "AxBxS", "AxC", "AxCxS", "BxC", "BxCxS", "AxBxC", "AxBxCxS", "T")
    result <- list(ANOVA.table=ANOVA.table, tbl=tbl)
    class(result) <- "SABC"
    return(result)
}
# ANOVA 表の print メソッド（SABC 関数が返すオブジェクト）
print.SABC <- function(obj) {
    x <- obj$ANOVA.table
    printf <- function(x, fmt) if (is.na(x)) "" else sprintf(fmt, x)
    x[,4] <- sapply(x[,4], printf, "%.5f")
    x[,5] <- sapply(x[,5], printf, "%.5f")
    print.data.frame(x)
}
# 推定周辺平均を出力する summary メソッド（SABC 関数が返すオブジェクト）
summary.SABC <- function(obj) {
    tbl <- obj$tbl
    d <- dim(tbl)
    n <- d[4]
    r <- d[1]
    q <- d[2]
    p <- d[3]
    means.a <- apply(tbl, c(4, 3), sum)/q/r
    means.b <- apply(tbl, c(4, 2), sum)/p/r
    means.c <- apply(tbl, c(4, 1), sum)/p/q
    mean <- c(colMeans(means.a), colMeans(means.b), colMeans(means.c))
    SE <- c(apply(means.a, 2, sd),
            apply(means.b, 2, sd),
            apply(means.c, 2, sd)) / sqrt(n)
    t95 <- qt(0.975, n - 1)
    EMM <- data.frame(Mean = mean, SE = SE,
                      LCL = mean - t95 * SE,
                      UCL = mean + t95 * SE)
    rownames(EMM) <- c(paste("a", 1:p, sep="="),
                       paste("b", 1:q, sep="="),
                       paste("c", 1:r, sep="="))
    print(EMM)
}
SG <- function(x)                                                       # データベクトル
{
        method <- "スミルノフ・グラブス検定"
        data.name <- paste(c("min(", "max("), deparse(substitute(x)), ") = ", range(x, na.rm=TRUE), sep="")
        x <- x[!is.na(x)]                                               # 欠損値を除く
        n <- length(x)                                                  # 標本サイズ
        t <- abs(range(x)-mean(x))/sd(x)                                # 最大のデータと最小のデータの両方について検定統計量を計算する
        p <- n*pt(sqrt((n-2)/((n-1)^2/t^2/n-1)), n-2, lower.tail=FALSE) # P 値も2通り計算される
        p <- sapply(p, function(p0) min(p0, 1))
        result <- list(method=method, parameter=c(df=n-2))
        result1 <- structure(c(result,  list(data.name=data.name[1], statistic=c(t=t[1]), p.value=p[1])), class="htest")
        result2 <- structure(c(result,  list(data.name=data.name[2], statistic=c(t=t[2]), p.value=p[2])), class="htest")
        return(structure(list(result1, result2), class="SG"))
}
# シェッフェの一対比較法 Scheffe's Paired Comparison
ScheffePairedComparison <- function(	A,			# 一対比較の結果を表す正方行列
					B,			# 得点ベクトル
					labels=NULL)		# 評価対象のラベル
{
	n <- (1+sqrt(1+8*nrow(A)))/2				# 評価対象の個数
	if (is.null(labels)) labels <- LETTERS[1:n]		# 評価対象名の補完
	AB <- A%*%B						# 従属変数ベクトルの作成
	x <- combn(n, 2)					# 独立変数行列の作成
	nc <- ncol(x)
	indep <- matrix(0, nc, n)
	indep[cbind(1:nc, x[1,])] <- 1
	indep[cbind(1:nc, x[2,])] <- -1
	ans <- lm(AB ~ indep[,2:n])				# 重回帰分析
	y <- coefficients(ans)					# 回帰係数がスコアになる
	names(y) <- labels
	return(structure(list(score=y, sorted.score=sort(y)), class="ScheffePairedComparison"))
}
# print メソッド
print.ScheffePairedComparison <- function(obj,			# ScheffePairedComparison が返すオブジェクト
				digits=5)			# 結果の表示桁数
{
	cat("\nスコア\n\n")
	print(round(obj$score, digits=digits))
	cat("\nソートされたスコア\n\n")
	print(round(obj$sorted.score, digits=digits))
}
# plot メソッド
plot.ScheffePairedComparison <- function(obj,			# ScheffePairedComparison が返すオブジェクト
					xlab="Score",		# 結果グラフの横軸名
					main="Scheffe's Paired Comparison",	# 結果グラフの表題
					file="")		# 結果グラフをファイル出力するときにファイル名
{
	if (file != "") pdf(file, width=540/72, height=160/72, onefile=FALSE)
	score <- obj$score
	plot(score, rep(0, length(score)), pch=19, xlab=xlab, main=main, , xaxt="n",
		xlim=range(pretty(score)), ylab="", yaxt="n", ylim=c(0,0.2),
		bty="n", xpd=TRUE)
	text(score, 0.0, names(score), pos=3)
	axis(1, pos=0)
	if (file != "") dev.off()	
}
# 連関比率法
SeasonalIndex <- function(	x,						# 4 半期ごとのデータベクトル
				xlab="", ylab="", main=NULL,			# 軸，図のラベル
				lx="bottomright", ly=NULL,			# legend の位置
				lty1=2, lty2=1,					# 線種
				pch1=1, pch2=19,				# マーカー
				label1="粗データ", label2="季節調整済みデータ")	# ラベル
{
	n <- length(x)
	y <- matrix(x/c(x[1], x[-n]), 4)
	mean1 <- rowMeans(y)
	mean2 <- mean1/exp(mean(log(mean1)))
	mean2[1] <- 1
	chain.index <- cumprod(mean2)						# 連鎖指数
	seasonal.index <- chain.index/mean(chain.index)				# 季節指数
	z <- as.vector(x/seasonal.index)					# 季節調整済みデータ
	plot(1:n, x, type="l", lty=lty1, xlab=xlab, ylab=ylab, main=main)
	points(1:n, x, pch=pch1)
	lines(1:n, z, lty=lty2)
	points(1:n, z, pch=pch2)
	legend(lx, ly, legend=c(label1, label2), lty=c(lty1, lty2), pch=c(pch1, pch2))
	return(list(seasonal.index=seasonal.index, corrected.data=z))
}
# シャーリー・ウィリアムズの方法による多重比較　
Shirley.Williams <- function(	data,				# データベクトル
				group,				# 群変数ベクトル
				method=c("up", "down"))		# 方法の指定
{
	OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
	data <- data[OK]
	group <- group[OK]
	method <- match.arg(method)				# 引数の補完
	func <- if (method == "down") min else max		# method により，後で使う関数を選ぶ
	ni <- table(group)					# 各群のデータ数
	a <- length(ni)						# 群の数
	s <- numeric(a-1)
	for (p in a:2) {
		select <- 1 <= group & group <= p		# 分析対象にするケースの選択・非選択ベクトル
		r <- rank(data[select])				# 順位付け
		g <- group[select]				# 分析対象とされたデータの群変数データ
		M <- func(cumsum(rev(tapply(r, g, sum))[-p])/cumsum(rev(ni[2:p])))
		N <- sum(ni[1:p])
		V <- (sum(r^2)-N*(N+1)^2/4)/(N-1)
		t <- (M-sum(r[group == 1])/ni[1])/sqrt(V*(1/ni[p]+1/ni[1]))
		s[p-1] <- ifelse(method == "down", -t, t)
	}
	t <- rev(s)
	names(t) <- a:2
	return(t)
}
# スティール・ドゥワス(Steel-Dwass)の方法による多重比較
Steel.Dwass <- function(data,						# データベクトル
			group)						# 群変数ベクトル
{
	OK <- complete.cases(data, group)				# 欠損値を持つケースを除く
	data <- data[OK]
	group <- group[OK]
	n.i <- table(group)						# 各群のデータ数
	ng <- length(n.i)						# 群の数
	t <- combn(ng, 2, function(ij) {
		i <- ij[1]
		j <- ij[2]
		r <- rank(c(data[group == i], data[group == j]))	# 群 i, j をまとめてランク付け
		R <- sum(r[1:n.i[i]])					# 検定統計量
		N <- n.i[i]+n.i[j]					# 二群のデータ数の合計
		E <- n.i[i]*(N+1)/2					# 検定統計量の期待値
		V <- n.i[i]*n.i[j]/(N*(N-1))*(sum(r^2)-N*(N+1)^2/4)	# 検定統計量の分散
		return(abs(R-E)/sqrt(V))				# t 値を返す
	})
	p <- ptukey(t*sqrt(2), ng, Inf, lower.tail=FALSE)		# P 値を計算
	result <- cbind(t, p)						# 結果をまとめる
	rownames(result) <- combn(ng, 2, paste, collapse=":")
	return(result)
}
# スティール（Steel）の方法による多重比較
Steel <- function(	data,					# データベクトル
			group)					# 群変数ベクトル
{
	get.rho <- function(ni)					# ρを計算する
	{
		k <- length(ni)
		rho <- outer(ni, ni, function(x, y) { sqrt(x/(x+ni[1])*y/(y+ni[1])) })
		diag(rho) <- 0
		sum(rho[-1, -1])/(k-2)/(k-1)
	}

	OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
	data <- data[OK]
	group <- factor(group[OK])				# 群変数データを factor に変関る
	ni <- table(group)					# 各群のデータ数
	a <- length(ni)						# 群の数
	control <- data[group == 1]				# 対照群のデータを取り出す
	n1 <- length(control)					# 対照群のデータ数
	t <- numeric(a)
	rho <- ifelse(sum(n1 == ni) == a, 0.5, get.rho(ni))	# ρを決める
	for (i in 2:a) {
		r <- rank(c(control, data[group == i]))		# 対照群と対照群をまとめて順位をつける
		R <- sum(r[1:n1])				# 検定統計量
		N <- n1+ni[i]					# 対照群と対照群のデータ数の合計
		E <- n1*(N+1)/2					# 検定統計量の期待値
		V <- n1*ni[i]/N/(N-1)*(sum(r^2)-N*(N+1)^2/4)	# 検定統計量の分散
		t[i] <- abs(R-E)/sqrt(V)			# 検定統計量（t 分布による漸近近似）
	}
	result <- cbind(t, rho)[-1,]				# 結果をまとめる
	rownames(result) <- paste(1, 2:a, sep=":")
	return(result)
}
# サーストンの一対比較法 Thurstone's Paired Comparison
ThurstonePairedComparison <- function(x)			# 一対比較の結果を表す正方行列
{
	nc <- ncol(x)						# 項目数
	stopifnot(nc == nrow(x))				# 正方行列でないと分析中止
	if (is.null(dimnames(x))) {				# 項目名がないときは補完する
		labels <- LETTERS[1:nc]
	}
	else if(is.null(colnames(x))) {
		labels <- rownames(x)
	}
	else {
		labels <- colnames(x)
	}
	n <- x+t(x)						# 対戦総数（引き分けとか試合数不足を考慮）
	diag(n) <- 1						# 0 による割り算が起きないように対角成分を調整
	x <- qnorm(x/n)						# 割合を求め，対応する Z スコアを求める
	diag(x) <- NA						# 対角は NA にする
	score <- rowMeans(x, na.rm=TRUE)			# 行和が求める答え
	names(score) <- labels					# スコアに項目名をつける
	return(structure(list(score=score, sorted.score=sort(score)), class="ThurstonePairedComparison"))
}
# print メソッド
print.ThurstonePairedComparison <- function(	obj,		# ThurstonePairedComparison が返すオブジェクト
						digits=5)
{
	cat("\nスコア\n\n")
	print(round(obj$score, digits=digits))
	cat("\nソートされたスコア\n\n")
	print(round(obj$sorted.score, digits=digits))
}
# plot メソッド
plot.ThurstonePairedComparison <- function(	obj,		# ThurstonePairedComparison が返すオブジェクト
						xlab="Score",	# 結果グラフの横軸名
						main="Thurstone's Paired Comparison",	# 結果グラフの表題
						file="")	# 結果グラフをファイル出力するときにファイル名
{
	if (file != "") pdf(file, width=540/72, height=160/72, onefile=FALSE)
	score <- obj$score
	plot(score, rep(0, length(score)), pch=19, xlab=xlab, main=main, xaxt="n",
		xlim=range(pretty(score)), ylab="", yaxt="n", ylim=c(0,0.2),
		bty="n", xpd=TRUE)
	text(score, 0.0, names(score), pos=3)
	axis(1, pos=0)
	if (file != "") dev.off()
}
# ライアンの方法とチューキーの方法による平均値の対比較
m.multi.comp <- function(	n,				# 標本サイズベクトル
				me,				# 平均値ベクトル
				s,				# 標準偏差ベクトル
				alpha=0.05,			# 有意水準
				method=c("ryan", "tukey"))	# 方法
{
	printf <- function(fmt, ...)
	{
		cat(sprintf(fmt, ...))
	}

	check <- function(s, b)					# 検定しようとしている二群が，それまでに有意でないとされた二群に挟まれているか
	{
		if (ns.n > 1) {
			for (i in 1:ns.n) {
				if (ns.s[i] <= s && s <= ns.b[i] && ns.s[i] <= b && b <= ns.b[i]) {
					return(FALSE) 		# 検定するまでもなく有意でないとする
				}
			}
		}
		return(TRUE)					# 検定しなくてはならない
	}

	k <- length(n)						# 群の数
	stopifnot(k == length(me), k == length(s), n > 0, floor(n) == n, s > 0)
	method <- match.arg(method)				# 引数の補完
	o <- order(me)						# 平均値の大きさの順位
	sn <- n[o]						# 並べ替えた標本サイズ
	sm <- me[o]						# 並べ替えた平均値
	ss <- s[o]						# 並べ替えた標準偏差
	nt <- sum(sn)						# 全体の標本サイズ
	mt <- sum(sn*sm)/nt					# 全体の平均値
	dfw <- nt-k						# 群内平方和の自由度
	vw <- sum((sn-1)*ss^2)/dfw				# 群内分散
	num.significant <- ns.n <- 0
	ns.s <- ns.b <- numeric(k*(k-1)/2)			# 有意でない群の記録用
	for (m in k:2) {					# 検定対象の選定
		for (small in 1:(k-m+1)) {
			big <- small+m-1
			if (check(small, big)) {
				t0 <- (sm[big]-sm[small])/sqrt(vw*(1/sn[big]+1/sn[small]))	# 検定統計量
				if (method == "ryan") { 					# Ryan の方法
					P <- pt(t0, dfw, lower.tail=FALSE)*2			# 有意確率
					nominal.alpha <- 2*alpha/(k*(m-1))			# 名義的有意水準
					result <- P <= nominal.alpha				# 検定結果
				}
				else { # Tukey の方法
					P <- ptukey(t0, m, dfw, lower.tail=FALSE)		# 有意確率
					result <- P <= alpha					# 検定結果
				}
				if (result) { 			# 有意であるとき
					num.significant <- 1
					printf("mean[%2i]=%7.5f vs. mean[%2i]=%7.5f : diff.= %7.5f, ",
						o[small], me[small], o[big], me[big], me[big]-me[small])
					if (method == "ryan") {
						printf("t=%7.5f : P=%7.5f, alpha'=%7.5f\n", t0, P, nominal.alpha)
					}
					else {
						printf("t=%7.5f : P=%7.5f\n", t0, P)
					}
				}
				else {				# 有意でないとき
					ns.n <- ns.n+1
					ns.s[ns.n] <- small
					ns.b[ns.n] <- big
				}
			}
		}
	}
	if (num.significant == 0) {				# 有意差のある群は一つもなかった
		print("Not significant at all.")
	}
}
# マン・ホイットニーの U 検定
U.test <- function(	x,				# 第一群の観測値ベクトルまたは，分割表データ（y=NULL)
			y = NULL,			# 第二群の観測値ベクトル
			correct = TRUE)			# 連続性の補正を行うかどうか
{
	method <- "マン・ホイットニーの U 検定"
	if (is.null(y)) {				# 2 × C 行列の分割表として与えられたとき
		if (nrow(x) != 2) stop("2 x C matrix is expected.")
		data.name <- paste(deparse(substitute(x)), "as 2 by C matrix")
		nc <- ncol(x)				# カテゴリー数
		y <- x[2,]				# 第二群の度数分布
		x <- x[1,]				# 第一群の度数分布
		tie <- x+y				# 合計した度数分布（同順位）
		n1 <- sum(x)				# 第一群のサンプルサイズ
		n2 <- sum(y)				# 第二群のサンプルサイズ
		n <- n1+n2				# 合計したサンプルサイズ
		rj <- c(0, cumsum(tie)[-nc])+(tie+1)/2	# カテゴリーに属するものの順位
		U1 <- n1*n2+n1*(n1+1)/2-sum(x*rj)	# 検定統計量
	}
	else {						# 2 つのデータベクトルとして与えられたとき
		data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
		x <- x[!is.na(x)]			# 欠損値を持つケースを除く
		y <- y[!is.na(y)]			# 欠損値を持つケースを除く
		n1 <- length(x)				# 第一群のサンプルサイズ
		n2 <- length(y)				# 第二群のサンプルサイズ
		n <- n1+n2				# 合計したサンプルサイズ
		xy <- c(x, y)				# 両群のデータを結合したもの
		r <- rank(xy)				# 順位のベクトル
		U1 <- n1*n2+n1*(n1+1)/2-sum(r[1:n1])	# 検定統計量
		tie <- table(r)				# 同順位の集計
	}
	U <- min(U1, n1*n2-U1)				# U 統計量
	V <- n1*n2*(n^3-n-sum(tie^3-tie))/12/(n^2-n)	# 同順位を考慮した分散
	E <- n1*n2/2					# 期待値
	Z <- (abs(U-E)-ifelse(correct, 0.5, 0))/sqrt(V)	# Z 値
	P <- pnorm(Z, lower.tail=FALSE)*2		# 両側 P 値
	return(structure(list(statistic=c(U=U, "E(U)"=E, "V(U)"=V, "Z-value"=Z), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# ウィリアムズの方法による多重比較を行う
Williams <- function(	data,				# データベクトル
			group,				# 群を表すベクトル
			method=c("up", "down"))		# 方法の指定
{
	OK <- complete.cases(data, group)		# 欠損値を持たないケースを選択する
	data <- data[OK]
	group <- group[OK]
	method <- match.arg(method)			# 引数の補完
	func <- if (method == "down") min else max	# 引数により，min か max かを選ぶ
	n.i <- tapply(data, group, length)		# 各群の例数
	sum.i <- tapply(data, group, sum)		# 群ごとの総和
	mean.i <- tapply(data, group, mean)		# 群ごとの平均値
	v.i <- tapply(data, group, var)			# 群ごとの不偏分散
	a <- length(n.i)				# 群の数	
	phi.e <- sum(n.i)-a				# 誤差分散の自由度
	v.e <- sum((n.i-1)*v.i)/phi.e			# 誤差分散
	t <- sapply(a:2,				# t 値を計算
		    function(p) (func(cumsum(rev(sum.i[2:p]))/cumsum(rev(n.i[2:p]))) - mean.i[1])/sqrt(v.e*(1/n.i[1]+1/n.i[p])))
	names(t) <- c(a:2)				# 名前を付ける
	return(list(phi.e=phi.e, t=if (method == "down") -t else t))
}
# 自己相関係数を計算する
# R にも用意されている
acf2 <- function(	x,	# 時系列データ
			k)	# ラグ
{
	n <- length(x)
	if (n < 3 || n-k < 2 || k < 1) {
		stop("invalid argument")
	}
	mean <- mean(x)
	num <- sum((x[1:(n-k)]-mean)*(x[(k+1):n]-mean))
	den <- var(x)*(n-1)
	return(num/den)
}
# 総当たり法による正準判別分析
all.candis <- function(	data,						# 説明変数のデータフレーム（データ行列）
			group)						# 群を表す変数（データフレームから指定するときには iris[,5] ではなく，iris[5] のように）
{
	BinConv <- function(nv)
	{
	       n <- 2^nv						# 独立変数を取り出す取り出し方
	        bincomb <- matrix(FALSE, nrow=n, ncol=nv)		# e1071 パッケージの bincombinations より
	        for (j in 1:nv) {
	                bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
	        }
	        bincomb <- bincomb[-1,]
		return(bincomb)
	}
	nv <- ncol(data)
	vname <- colnames(data)						# 変数名（なければ作る）
	if (is.null(vname)) {
		vname <- colnames(data) <- paste("x", 1:nv, sep="")
	}
	gname <- names(group)
	if (is.null(gname)) {						# group を，データフレームから iris[, 5] のようにすると，名前がなくなる
		gname <- ""						# iris[5] のように 1 列のデータフレームとして指定すること
	}
	group <- factor(as.matrix(group))
	ok <- complete.cases(data, group)				# 欠損値を含まないケースだけを対象にする
	data <- data[ok,]
	group <- group[ok]
	n <- nrow(data)
	bincomb <- BinConv(nv)
	nr <- nrow(bincomb)
	correct <- numeric(nr)						# 正準判別関数の正判別率を基準とする
	for (i in 1:nr) {
		dat <- data[, bincomb[i,], drop=FALSE]
		a <- candis(dat, group)					# 別途用意してある candis.R で定義
		correct[i] <- sum(a$classification == group)/n		# 正判別率を記録しておく
	}
	ans <- data.frame(correct, bincomb)
	colnames(ans) <- c("correct rate", vname)
	return(structure(list(ans=ans, name=vname, gname=gname), class="all.candis"))
}
# print メソッド
print.all.candis <- function(obj)
{
	ans <- obj$ans
	name <- obj$name
	gname <- obj$gname
	o <- order(ans[, 1], decreasing=TRUE)
	ans <- ans[o,]
	nc <- ncol(ans)
	cat("\ncorrect rate   formula\n")
	for (i in 1:nrow(ans)) {
		cat(sprintf("%10.5f     %s ~ %s\n",ans[i, 1], gname, paste(name[as.matrix(ans[i, 2:nc])], collapse=" + ")))
	}
	invisible(ans)							# 結果をソートしただけのものを返す
}
# 総当たり法による線形判別分析
all.disc <- function(	data,						# 説明変数のデータフレーム（データ行列）
			group)						# 群を表す変数（データフレームから指定するときには iris[,5] ではなく，iris[5] のように）
{
	BinConv <- function(nv)
	{
	       n <- 2^nv						# 独立変数を取り出す取り出し方
	        bincomb <- matrix(FALSE, nrow=n, ncol=nv)		# e1071 パッケージの bincombinations より
	        for (j in 1:nv) {
	                bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
	        }
	        bincomb <- bincomb[-1,]
		return(bincomb)
	}
	nv <- ncol(data)
	vname <- colnames(data)						# 変数名（なければ作る）
	if (is.null(vname)) {
		vname <- colnames(data) <- paste("x", 1:nv, sep="")
	}
	gname <- names(group)
	if (is.null(gname)) {						# group を，データフレームから iris[, 5] のようにすると，名前がなくなる
		gname <- ""						# iris[5] のように 1 列のデータフレームとして指定すること
	}
	group <- factor(as.matrix(group))
	ok <- complete.cases(data, group)				# 欠損値を含まないケースだけを対象にする
	data <- data[ok,]
	group <- group[ok]
	n <- nrow(data)
	bincomb <- BinConv(nv)
	nr <- nrow(bincomb)
	correct <- numeric(nr)						# 正準判別関数の正判別率を基準とする
	for (i in 1:nr) {
		dat <- data[, bincomb[i,], drop=FALSE]
		a <- disc(dat, group)					# 別途用意してある disc.R で定義
		correct[i] <- sum(a$correct)/n		# 正判別率を記録しておく
	}
	ans <- data.frame(correct, bincomb)
	colnames(ans) <- c("correct rate", vname)
	return(structure(list(ans=ans, name=vname, gname=gname), class="all.disc"))
}
# print メソッド
print.all.disc <- function(obj)
{
	ans <- obj$ans
	name <- obj$name
	gname <- obj$gname
	o <- order(ans[, 1], decreasing=TRUE)
	ans <- ans[o,]
	nc <- ncol(ans)
	cat("\ncorrect rate   formula\n")
	for (i in 1:nrow(ans)) {
		cat(sprintf("%10.5f     %s ~ %s\n",ans[i, 1], gname, paste(name[as.matrix(ans[i, 2:nc])], collapse=" + ")))
	}
	invisible(ans)							# 結果をソートしただけのものを返す
}
# 総当たり法によるロジスティック回帰を行う
# 　データフレームには，分析に使用する独立変数と従属変数のみを含むこと。
# 　また，従属変数は最終列に置くこと。
#
all.logistic <- function(	df,				# データフレーム（独立変数）
				limit=10,			# 独立変数の個数の上限（数が多いと計算時間が指数的に増える）
				max.p=NULL,			# モデルに組み入れる独立変数の数の上限
				min.p=max.p)			# モデルに組み入れる独立変数の数の下限
{
	df <- subset(df, complete.cases(df))			# 欠損値を持つケースを除く
	nv <- ncol(df)-1					# 独立変数の個数
	if (is.null(max.p)) {					# 完全な総当たり
		if (nv > limit) {				# limit より多いと分析を中断する
			stop(paste("独立変数が", limit,
			"個以上である（多すぎる）。\n",
			"limit 引数で変更できる", paste=""))
		}
		n <- 2^nv					# 独立変数を取り出す取り出し方
		bincomb <- matrix(FALSE, nrow=n, ncol=nv)	# e1071 パッケージの bincombinations より
		for (j in 1:nv) {
			bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
		}
		bincomb <- bincomb[-1,]
	}
	else {							# モデルに組み入れる独立変数の数を制限する
		if (max.p < 1 || max.p > nv || min.p < 1 || min.p > nv || min.p > max.p) {
			stop("モデルに取り込む独立変数の数の指定が不正")
		}
		bincomb <- NULL
		for (np in min.p:max.p) {
			bc <- combn(nv, np, function(x) 1:nv %in% x)
			bincomb <- rbind(bincomb, t(bc))
		}
	}
	n <- nrow(bincomb)
	name <- names(df)					# 変数名を取り出す
	depname <- name[nv+1]
	name <- name[1:nv]
	result3 <- character(n)					# 数値型ベクトル確保
	result1 <- result2 <- numeric(n)			# 数値型ベクトル確保
	for (i in 1:n) {					# 独立変数の全ての組み合わせについて，
		str <- name[bincomb[i,]]			# どの独立変数が使われるかを割り出す
		form <- reformulate(str, depname)		# モデル式を作る（"formula" クラス）
		ans <- glm(form, df, family=binomial)		# ロジスティック回帰分析の結果
		result <- summary(ans)
		result1[i] <- result$deviance			# deviance
		result2[i] <- result$aic			# AIC
		temp <- as.character(form)			# モデル式を文字列に変換
		result3[i] <- paste(temp[2], "~", temp[3])	# モデル式を記録
	}
	return(structure(list(deviance=result1, AIC=result2, formula=result3),
	class="all.logistic"))
}
# print メソッド
print.all.logistic <- function( obj,				# "all.logistic" クラスのオブジェクトをプリント
				sort.by=c("deviance", "AIC"),	# 結果を何で並べ替えるかを指示
				models=20) 			# 良い方から何番目まで出力するか
{
	result <- data.frame(obj$deviance, obj$AIC, obj$formula)
	sort.by <- match.arg(sort.by)
	o <- order(switch (sort.by, "deviance"=result[,1], "AIC"=result[,2]))
	result <- result[o,]
	cat(sprintf("\n%12s %12s   %s\n", "deviance", "AIC", "Formula"))	# 表頭
	models <- min(models, nrow(result))
	apply(result[1:models,], 1, function(x)			# 各行の出力
		cat(sprintf("%12.5f %12.5f   %s\n",
		as.double(x[1]), as.double(x[2]), x[3])))
	invisible(result)
}
# 総当たり法による二次の判別分析
all.quad.disc <- function(	data,						# 説明変数のデータフレーム（データ行列）
			group)						# 群を表す変数（データフレームから指定するときには iris[,5] ではなく，iris[5] のように）
{
	BinConv <- function(nv)
	{
	       n <- 2^nv						# 独立変数を取り出す取り出し方
	        bincomb <- matrix(FALSE, nrow=n, ncol=nv)		# e1071 パッケージの bincombinations より
	        for (j in 1:nv) {
	                bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
	        }
	        bincomb <- bincomb[-1,]
		return(bincomb)
	}
	nv <- ncol(data)
	vname <- colnames(data)						# 変数名（なければ作る）
	if (is.null(vname)) {
		vname <- colnames(data) <- paste("x", 1:nv, sep="")
	}
	gname <- names(group)
	if (is.null(gname)) {						# group を，データフレームから iris[, 5] のようにすると，名前がなくなる
		gname <- ""						# iris[5] のように 1 列のデータフレームとして指定すること
	}
	group <- factor(as.matrix(group))
	ok <- complete.cases(data, group)				# 欠損値を含まないケースだけを対象にする
	data <- data[ok,]
	group <- group[ok]
	n <- nrow(data)
	bincomb <- BinConv(nv)
	nr <- nrow(bincomb)
	correct <- numeric(nr)						# 正準判別関数の正判別率を基準とする
	for (i in 1:nr) {
		dat <- data[, bincomb[i,], drop=FALSE]
		a <- quad.disc(dat, group)					# 別途用意してある quad.disc.R で定義
		correct[i] <- sum(a$correct)/n		# 正判別率を記録しておく
	}
	ans <- data.frame(correct, bincomb)
	colnames(ans) <- c("correct rate", vname)
	return(structure(list(ans=ans, name=vname, gname=gname), class="all.quad.disc"))
}
# print メソッド
print.all.quad.disc <- function(obj)
{
	ans <- obj$ans
	name <- obj$name
	gname <- obj$gname
	o <- order(ans[, 1], decreasing=TRUE)
	ans <- ans[o,]
	nc <- ncol(ans)
	cat("\ncorrect rate   formula\n")
	for (i in 1:nrow(ans)) {
		cat(sprintf("%10.5f     %s ~ %s\n",ans[i, 1], gname, paste(name[as.matrix(ans[i, 2:nc])], collapse=" + ")))
	}
	invisible(ans)							# 結果をソートしただけのものを返す
}
# クロンバックのα信頼性係数
alpha <- function(	x,			# 必要な変数のみからなるデータ行列（データフレームまたは行列）
			detail=FALSE)		# 詳細情報を計算するかどうか
{
	alpha0 <- function(x)
	{
		k <- ncol(x)			# 変数の個数
		VarCovMat <- var(x)		# 分散共分散行列
		Sy2 <- sum(VarCovMat)		# 合計点の不偏分散（var(rowSums(x)) と同じ）
		Sj2 <- sum(diag(VarCovMat))	# 各変数の不偏分散の和
		return(k/(k-1)*(1-Sj2/Sy2))
	}
	x <- as.matrix(x)			# 行列にする
	x <- subset(x, complete.cases(x))	# 欠損値を持つケースを除く
	k <- ncol(x)				# 変数の個数（列数）
       if (is.null(colnames(x))) {		# 変数名が無いときには仮の名前を付ける
                colnames(x) <- paste("X", 1:k, sep="")
        }
        vnames <- colnames(x)			# 変数名を記録
	stopifnot(k > 1,			# 列数は2列以上であること
		  nrow(x) > 1)			# 行数も2行以上であること…さもなくば，中止
	alpha <- alpha0(x)
	alpha2 <- cor2 <- R2 <- numeric(k)
	if (detail == TRUE && k >= 3) {		# print メソッドの注を参照
		z <- rowSums(x)
		for (i in 1:k) {
			x2 <- x[, -i]
			alpha2[i] <- alpha0(x2)
			cor2[i] <- cor(x[, i], rowSums(x2))
			R2[i] <- summary(lm(x[,i] ~ x2))$r.squared
		}
		result2 <- data.frame(alpha_=alpha2, r_=cor2, R_=R2)
		rownames(result2) <- vnames
		return(structure(list(alpha=c(alpha=alpha), result2=result2), class="alpha"))
	}
	else {
		class(alpha) <- "alpha"
		return(alpha=alpha)
	}
}
# print メソッド
print.alpha <- function(obj, digits=5)
{
	fmt <- sprintf("alpha = %%.%if\n", digits)
	if (length(obj) == 1) {
		cat(sprintf(fmt, obj))
	}
	else {
		cat(sprintf(fmt, obj$alpha), "\n")
		print(obj$result2, digits=digits)
		cat("alpha_ : それぞれの変数を除いたときの alpha\n")
		cat("r_     : それぞれの変数とその変数を除いたときの合計値との相関係数\n")
		cat("R_     : それぞれの変数をその変数以外の変数で予測したときの決定係数\n")
	}
}
# 多角形の面積を計算する
Area <- function(xy)					# 座標データ行列（1列目がx座標，2列目がy座標）
{							# 時計回りの順番で用意する
	n <- nrow(xy)
	x <- xy[,1]
	# 0.5 * Σ ((x[2], x[3], ..., x[n], x[1]) - (x[n], x[1], x[2], ..., x[n-1])) * xy[,2]
	# sum(sapply(1:n, function(i) x[i%%n+1]-x[(i+n-2)%%n+1])*xy[,2])/2
	sum((c(x[-1], x[1])-c(x[n], x[-n]))*xy[,2])/2
}
# 度数分布表から基礎統計量を求める
basic.stat <- function(	x,								# 級限界のベクトル
			f)								# 度数のベクトル
{
	w <- diff(x[1:2])								# 区間幅
	stopifnot(all(diff(x) == w))							# 区間は等間隔でなければならない
	stopifnot(length(x) == length(f))						# ベクトルの長さは同じでなければならない
	x <- x+w/2									# 級中心ベクトルに変換
	n <- sum(f)									# サンプルサイズ
	m <- sum(f*x)/n									# 平均値
	v <- sum(f*(x-m)^2)/(n-1)							# 不偏分散
	SD <- sqrt(v)									# 標準偏差
	CV <- SD/m*100									# 変動係数
	g1 <- n*sum(f*(x-m)^3)/(n-1)/(n-2)/SD^3						# 歪度（不偏推定値）
	g2 <- n*(n+1)*sum(f*(x-m)^4)/(n-1)/(n-2)/(n-3)/SD^4-3*(n-1)^2/(n-2)/(n-3)	# 尖度
	result <- list(n=n, mean=m, variance=v, sd=SD, g1=g1, g2=g2, CV=CV)		# リストで返す
	class(result) <- c("basic.stat", "list")
	return(result)
}
print.basic.stat <- function(x)								# basic.stat のプリント・メソッド
{
	cat("標本の大きさ =", x$n); cat("\n")
	cat("算術平均値　 =", x$mean); cat("\n")
	cat("不偏分散　　 =", x$variance); cat("\n")
	cat("標準偏差　　 =", x$sd); cat("\n")
	cat("歪度　　　　 =", x$g1); cat("\n")
	cat("尖度　　　　 =", x$g2); cat("\n")
	cat("変動係数　　 =", x$CV); cat("\n")
}
# 二項分布への適合度の検定
binomdist <- function(	d,					# 度数分布ベクトル
			x,					# 階級値ベクトル
			size)					# 試行回数
{
	data.name <- paste( "\n度数分布ベクトル：", deparse(substitute(d)),
			"\n階級値ベクトル：　", deparse(substitute(x)),
			"\n試行回数：　　　　", deparse(substitute(size)), sep="")
	method <- "二項分布への適合度の検定"
	k <- length(d)						# 階級数
	if (k != length(x)) {
		stop("度数ベクトル d と階級値ベクトル x の長さが違います")
	}
	else if (any(floor(d) != d)) {
		stop("度数ベクトル中に小数値があります")
	}
	else if (any(d < 0)) {
		stop("度数ベクトル中に負の値があります")
	}
	else if (any(x > size)) {
		stop("階級値ベクトル中に試行数より大きい数値があります")
	}
	else if (any(x < 0)) {
		stop("階級値ベクトル中に負の数値があります")
	}
	else if (any(floor(x) != x)) {
		stop("階級値ベクトル中に小数値があります")
	}
	n <- sum(d)						# サンプルサイズ
	prob <- sum(d*x)/n/size					# 母比率
	p <- dbinom(x, size, prob)				# 二項分布の確率
	e <- n*p						# 期待値
	table <- data.frame(x, d, p, e)
	rownames(table) <- paste("c-", x, sep="")
	while (e[1] < 1) {					# 1 未満のカテゴリーの併合
		d[2] <- d[2]+d[1]
		e[2] <- e[2]+e[1]
		d <- d[-1]
		e <- e[-1]
		k <- k-1
	}
	while (e[k] < 1) {					# 1 未満のカテゴリーの併合
		d[k-1] <- d[k-1]+d[k]
		e[k-1] <- e[k-1]+e[k]
		d <- d[-k]
		e <- e[-k]
		k <- k-1
	}
	chisq <- sum((d-e)^2/e)					# カイ二乗統計量
	df <- k-2						# 自由度
	p <- pchisq(chisq, df, lower.tail=FALSE)		# P 値
	names(chisq) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=chisq, parameter=df,
		p.value=p, estimate=c(n=n, probability=prob), method=method,
		data.name=data.name, table=table),
		class=c("htest", "binomdist")))
}
# summary メソッド（適合度に関する結果を表示する）
summary.binomdist <- function(	obj,				# binomdist が返すオブジェクト
				digits=5)			# 表示桁数
{
	table <- obj$table
	colnames(table) <- c("級", "度数", "確率", "期待値")
	cat("\n適合度\n\n")
	print(table, digits=digits, row.names=FALSE)
}
# plot メソッド（観察度数と理論度数を図示する）
plot.binomdist <- function(	obj,				# binomdist が返すオブジェクト
				...)				# barplot 関数へ渡す引数
{
	table <- obj$table
	nr <- nrow(table)
	barplot(table$d, space=0, ...)				# 観察度数を barplot で描く
	old <- par(xpd=TRUE)
	points(1:nr-0.5, table$e, pch=3)			# は理論度数を，記号 + で示す
	text(1:nr-0.5, -strheight("H"), table$x)		# 階級表示
	par(old)
}
# 分法により一変数間数 f(x)=0 の解を求める
bisection <- function(	func,					# 関数　たとえば function(x) (x+6.7)*(x-3.4) のようなもの
			lower, upper,				# 解を求める区間
			ndiv=50,				# 区間を細区分する個数
			epsilon=1e-14,				# 解の許容誤差
			max.iteration=500)			# 反復回数の上限
{
	printf <- function(fmt, ...) cat(sprintf(fmt, ...))	# 書式付きの print
	bisec2 <- function(func, lower, upper)			# 区間の再区分について解を探索する下請け関数
	{
		yl <- func(lower)				# 区間の左端における関数値
		yu <- func(upper)				# 区間の右端における関数値
		if (yl*yu > 0) {				# 同じ符号のとき，この区間には解がない（区間の設定不良）
			return(1)				# 戻り値として 1 を返す（利用しなくてもよい）
		}
		for (i in 1:max.iteration) {			# 繰り返し上限まで収束計算を続ける
			mid <- (lower+upper)/2			# そのときの区間の中点の値
			ym <- func(mid)				# 中点における関数値
			if (abs(ym) < epsilon) {		# たまたま解になったら結果を書き出す
				printf("ans=%g\n", mid)
				return(0)
			}
			else if (yu*ym > 0) {			# 区間の右端での関数値と符号が同じなら，
				upper <- mid			# 中点を区間の右端にする
				yu <- ym			# 新しい区間の右端での関数値を設定する
			}
			else {					# 区間の左端での関数値と符号が同じなら，
				lower <- mid			# 中点を区間の左端にする
				yl <- ym			# 新しい区間の左端での関数値を設定する
			}
			if (abs(upper-lower) < epsilon) {	# 区間幅が誤差範囲内になったら中点を解とする
				printf("ans=%g\n", (lower+upper)/2)
				return(0)
			}
		}
	}

	x <- seq(lower, upper, length.out=ndiv)			# 区間を細区分して，
	for (i in 1:(ndiv-1)) {					# それぞれの細区分に対して，
		bisec2(func, x[i], x[i+1])			# 二分法により解を求める下請け関数を呼ぶ
	}
}
# Excel の二変量統計関数
correl <- cor						# 単に名前の違い。ただし，cor(x, y) の場合のみ

covar <- function(x, y)					# 共分散（不偏共分散ではない）
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	n <- length(x)					# ケース数
	var(x, y)*(n-1)/n				# Excel は共変動をデータ数（n）で割ったものとして定義している
}

forecast <- function(data, y, x)			# R の特徴で，data はスカラーでもベクトルでもかまわない
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	mean(y)+var(x, y)/var(x)*(data-mean(x))		# y，x について回帰直線を求め，独立変数の値が data のときの予測値を求める
}

growth <- function(y, x, data, one = FALSE)		# R の特徴で，data はスカラーでもベクトルでもかまわない
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	stopifnot(all(y > 0))				# 従属変数はすべて正の値を取らなければならない
	y <- log(y)					# y = a*b^x の，両辺の対数をとり，直線回帰に持ち込む
	if (one) {					# y，x について y = b^x にあてはめ，独立変数の値が data のときの予測値を求める
		b <- sum(x*y)/sum(x^2)			# 原点を通る直線
		const <- 0
	}
	else {						# y，x について y = a*b^x にあてはめ，独立変数の値が data のときの予測値を求める
		b <- var(x, y)/var(x)			# 切片を持つ直線
		const <- mean(y)-b*mean(x)
	}
	exp(const+b*data)				# 予測値を計算し，その指数をとって元の尺度に戻す
}

intercept <- function(y, x)
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	mean(y)-var(x, y)/var(x)*mean(x)		# 回帰直線の切片
}

logest <- function(y, x, one = FALSE)
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	stopifnot(all(y > 0))				# 従属変数はすべて正の値を取らなければならない
	y <- log(y)					# y = a*b^x の，両辺の対数をとり，直線回帰に持ち込む
	if (one) {					# y，x について y = b^x にあてはめ，b を求める
		b <- sum(x*y)/sum(x^2)			# 原点を通る直線
		const <- 0
	}
	else {						# y，x について y = a*b^x にあてはめ，a と bを求める
		b <- var(x, y)/var(x)			# 切片を持つ直線
		const <- mean(y)-b*mean(x)
	}
	list(model="a*b^x", result=c(a=exp(const), b=exp(b)))
}

pearson <- cor						# 単に名前の違い

rsq <- function(y, x)
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	cor(y, x)^2					# 単相関係数の二乗
}

slope <- function(y, x, zero = FALSE)
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	ifelse(zero, sum(x*y)/sum(x^2),			# zero=TRUE  のときには，原点を通る回帰直線の傾きを計算する
		var(x, y)/var(x))			# zero=FALSE のときには，切片を持つ回帰直線の傾きを計算する
}

steyx <- function(y, x)
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	n <- length(x)
	sqrt((var(y)-var(x, y)^2/var(x))*(n-1)/(n-2))	# 回帰直線において，回帰からの標本標準偏差（残差に対する平均平方の平方根）を計算する
}

trend <- function(y, x, data, zero = FALSE)		# R の特徴で，data はスカラーでもベクトルでもかまわない
{
	OK <- complete.cases(x, y)			# 二変数ともに欠損値を持たないケースを選択する
	x <- x[OK]
	y <- y[OK]
	
	ifelse(zero, sum(x*y)/sum(x^2)*data,		# zero=TRUE  のときには，原点を通る回帰直線による予測値を計算する
		intercept(y, x)+var(x, y)/var(x)*data)	# zero=FALSE のときには，切片を持つ回帰直線による予測値を計算する
}
# 母平均の検定・推定
boheikin <- function(	n,				# 標本サイズ
			xbar,			# 標本平均
			U=NULL,			# 標本不偏分散（母分散が未知の場合に指定する）
			mu=0,			# 母平均（信頼区間だけを求めるときには不要）
			sigma2=NULL,		# 母分散（母分散が既知の場合に指定する）
			conf.level=0.95)	# 信頼区間の信頼率
{
	if (!is.null(U)) {			# 母分散が未知のとき
		data.name <- sprintf("n = %s, mean = %s, variance = %s, μ = %s", n, xbar, U, mu)
		method <- "二次データによる母平均の検定と推定（母分散が未知のとき）"
		t <- abs(xbar-mu)/sqrt(U/n)
		df <- n-1
		p <- pt(t, df, lower.tail=FALSE)*2
		q <- qt(0.5-conf.level/2, df)
		conf.int <- xbar+c(q, -q)*sqrt(U/n)
		attr(conf.int, "conf.level") <- conf.level
		names(t) <- "t"
		names(df) <- "df"
		return(structure(list(statistic=t, parameter=df, p.value=p,
			conf.int=conf.int, method=method, data.name=data.name),
			class="htest"))
	}
	else if (!is.null(sigma2)) {		# 母分散が既知のとき
		data.name <- sprintf("n = %s, mean = %s, μ = %s, σ2 = %s", n, xbar, mu, sigma2)
		method <- "二次データによる母平均の検定と推定（母分散が既知のとき）"
		z <- abs(xbar-mu)/sqrt(sigma2/n)
		p <- pnorm(z, lower.tail=FALSE)*2
		q <- qnorm(0.5-conf.level/2)
		conf.int <- xbar+c(q, -q)*sqrt(sigma2/n)
		attr(conf.int, "conf.level") <- conf.level
		names(z) <- "Z-value"
		return(structure(list(statistic=z, p.value=p,
			conf.int=conf.int, method=method, data.name=data.name),
			class="htest"))
	}
	else {
		stop("不偏分散か母分散かどちらも NULL では計算できません")
	}
}
# 母平均の推定
boheikin.cl <- boheikin
# 母平均の検定・推定
boheikin.test <- boheikin
#####
#
# 複数のカテゴリー変数により多元分類を行い，各群の平均値，標準偏差を求め，必要なら平均値・代表値の差の検定を行う
#
#####

breakdown <- function(i,							# 分析対象の変数が入っているデータフレーム上の列番号または変数名ベクトル
		      j,							# 群を表す変数が入っているデータフレーム上の列番号または変数名ベクトル
		      df,							# データフレーム
		      latex=TRUE,						# LaTeX 形式で集計表を出力する（デフォルトは LaTeX 形式）
                      captions=NULL,						# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
                      labels=NULL,						# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
		      test=c("none", "parametric", "non-parametric"),		# デフォルト none では検定を行わない。検定を行うときはその種類を指定する
		      statistics=c("mean", "median"),				# （平均値，標準偏差）を求めるか（中央値，四分偏差）を求めるかを指定する
		      var.equal=FALSE,						# t-test, oneway の場合に等分散性を仮定するかどうかの引数
		      digits=3,							# 平均値，標準偏差を表示するときの小数点以下の桁数
                      output="",						# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
		      encoding=getOption("encoding"))				# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

	SIQ <- function(x) return(diff(fivenum(x)[c(2,4)]))                     # 四分偏差を求める下請け関数

	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

	test <- match.arg(test)							# 引数が省略形で与えられたときに，正確な名前をとる
	statistics <- match.arg(statistics)					# 引数が省略形で与えられたときに，正確な名前をとる
	if (statistics == "mean") {
		MEAN <- mean							# 位置の母数を求める関数：平均値
		SD <- sd							# 散布度を求める関数：標準偏差
		M.str <- "平均値"
		S.str <- "標準偏差"
	}
	else {
		MEAN <- median							# 位置の母数を求める関数：中央値
		SD <-  SIQ							# 散布度を求める関数：四分偏差
		M.str <- "中央値"
		S.str <- "四分偏差"
	}
	format <- paste("%.", digits, "f", sep="")				# 小数点以下 3 桁で出力する書式
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}
	index <- 0
	for (k in i) {								# ベクトルで指定されたすべての変数について分析する
		ok <- complete.cases(df[,k], df[,j])				# 欠損値を持たないケースを特定
		df2 <- df[ok,]							# 欠損値を持つケースを除く
		nl <- length(j)							# 何元分類にwなるか
		lst <- if (nl == 1) list(df2[,j])
		       else as.list(append(NULL, df2[,j]))
		x <- df2[, k]							# 分析対象変数
		nt <- length(x)							# 全体のデータ数
		mt <- MEAN(x)							# 全体の平均値
		st <- SD(x)							# 全体の標準偏差
		nms <- as.matrix(cbind(aggregate(df2[,k], lst, length),		# 各群のデータ数
				       aggregate(df2[,k], lst, MEAN)[,nl+1],	# 各群の平均値
				       aggregate(df2[,k], lst, SD)[,nl+1]))	# 各群の標準偏差を取り出し，行列形式にする
		nr <- nrow(nms)							# 行数（分類数）
		nc <- ncol(nms)							# 列数（分類詳細と，データ数，平均値，標準偏差）
		str <- paste(colnames(df2)[j], collapse="，")			# 多元分類に使われる変数名のリスト
		if (latex) {							# LaTeX 形式で集計結果を出力する
			index <- index+1
			cat("\n\\begin{table}[htbp]\n", file=output)		# \begin{table}[htbp]
			if (is.null(captions)) {
				cat(sprintf("\\caption{%s別の%sの集計}\n",	# \caption{○○別の□□の集計}
				    str, colnames(df2)[k]), file=output)
			}
			else {
				cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
			}
			if (!is.null(labels)) {
				cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
			}
			cat("\\centering\n", file=output)			# \centering
			cat("\\begin{tabular}{", rep("l", nc-3),		# \begin{tabular}{l…ccc} \hline
			    "ccc} \\hline\n", sep="", file=output)
			cat(rep("&", nc-3),					# 分類に使用する変数分の &
			    sprintf(" \\multicolumn{3}{c}{%s}\\\\ \\cline{%i-%i}\n", # 最後の3列を使って分析対象の変数名と罫線
			    colnames(df2)[k], nc-2, nc), file=output)
			cat(colnames(df2)[j], "データ数", M.str, S.str,		# 分類変数名 & データ数 & 平均値 & 標準偏差
			    sep=" & ", file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			for (l in 1:nr) {					# 各分類ごとに，
				cat(nms[l,1:(nc-2)],				# 分類基準，データ数，平均値，標準偏差
				    sprintf(format, as.numeric(nms[l, (nc-1):nc])), sep=" & ", file=output)
				cat("\\\\", file=output)			# \\
				if (l == nr) cat("\\hline\n", file=output)
				else cat("\n", file=output)			# 最後の行なら \hline
			}
			cat(sprintf("\\multicolumn{%i}{l}{全体}", nc-3),	# 全体について，データ数，平均値，標準偏差
			    nt, sprintf(format, mt), sprintf(format, st), sep=" & ", file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			cat("\\end{tabular}\n", file=output)			# \end{tabular}
		}
		else {							# tab で区切って出力する
			cat("\n表　", str, "別の", colnames(df2)[k], "の集計",	# 表　○○別の□□の集計
			    sep="", file=output)
			cat("\n", colnames(df2)[k], sep="\t", file=output,	# 分析対象変数
			    fill=TRUE)
			cat(colnames(df2)[j], "データ数", M.str, S.str,		# 分類変数名　データ数　平均値　標準偏差
			    sep="\t", file=output, fill=TRUE)
			for (l in 1:nr) {					# 各分類ごとに，
				cat(nms[l,1:(nc-2)], sprintf(format,		# 分類基準，データ数，平均値，標準偏差
				    as.numeric(nms[l, (nc-1):nc])), sep="\t",
				    file=output, fill=TRUE)
			}
			cat("全体", nt, sprintf(format, mt),			# 全体について，データ数，平均値，標準偏差
			    sprintf(format, st), sep="\t", file=output, fill=TRUE)
		}
		if (nr >= 2 && nl == 1) {					# 一元分類のときのみ検定を行う
			if (latex && test != "none") {
				cat("\\\\ \\noindent\n", file=output)
			}
			if (nr == 2 && test == "parametric") {			# 2 群の場合には t.test を呼ぶ
				res <- t.test(x~df2[,j], var.equal=var.equal)
				cat(sprintf(if (latex) "$t$値 = %.3f, 自由度 = %.3f, $P$値 = %.3f\n"
					    else "t 値 = %.3f, 自由度 = %.3f, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
			else if (nr >= 3 && test == "parametric") {		# 3 群以上の場合には oneway.test を呼ぶ
				res <- oneway.test(x~df2[,j], var.equal=var.equal)
				cat(sprintf(if (latex) "$F$値 = %.3f, 自由度 = (%i, %.3f), $P$値 = %.3f\n"
					    else "F 値 = %.3f, 自由度 = (%i, %.3f), P 値 = %.3f\n",
					    res$statistic, res$parameter[1], res$parameter[2], res$p.value), file=output)
			}
			else if (nr == 2 && test == "non-parametric") {		# 2 群以上の場合には wilcox.test を呼ぶ
				res <- wilcox.test(x~df2[,j])			# マン・ホイットニーの U 検定
				cat(sprintf(if (latex) "$U$ = %.3f, $P$値 = %.3f\n"
					    else "U = %.3f, P 値 = %.3f\n",
					    res$statistic, res$p.value), file=output)
			}
			else if (nr >= 3 && test == "non-parametric") {		# 3 群以上の場合には kruskal.test を呼ぶ
				res <- kruskal.test(x~df2[,j])
				cat(sprintf(if (latex) "$\\chi^2_{kw}$値 = %.3f, 自由度 = %i, $P$値 = %.3f\n"
					    else "カイ二乗値(kw) = %.3f, 自由度 = %i, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
		}
		if (latex) {							# LaTeX 形式で集計結果を出力する場合は，
			cat("\\end{table}\n", file=output)			# \end{table}
		}
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
# 正準判別分析
candis <- function(	data,							# 説明変数データ行列
			group)							# 群変数
{
	vnames <- colnames(data)
	group <- as.factor(as.matrix(group))					# 群を factor に変換
	OK <- complete.cases(data, group)					# 欠損値を持つケースを除く
	data <- as.matrix(data[OK,])
	colnames(data) <- vnames
	group <- group[OK]
	p <- ncol(data)								# 変数の個数
	n <- nrow(data)								# ケース数
	n.i <- table(group)							# 各群のケース数
	k <- length(n.i)							# 群の個数
	group.means <- matrix(unlist(by(data, group, colMeans)), p)		# 各変数の群ごとの平均値
	grand.means <- colMeans(data)						# 各変数の全体の平均値
	means <- cbind(grand.means, group.means)				# 結果表示のためにまとめる
	F.value <- apply(data, 2, function(x) oneway.test(x ~ group, var.equal=TRUE)$statistic)
	df1 <- k-1
	df2 <- n-k
	p.value <- format.pval(pf(F.value, df1, df2, lower.tail=FALSE))
	wilksFromF <- df2/(F.value*df1+df2)
	univariate <- data.frame(wilksFromF, F.value, rep(df1, p), rep(df2, p), p.value)
	ss <- split(data.frame(data), group)
	w <- Reduce("+", lapply(ss, function(x) var(x)*(nrow(x)-1)))		# 群内変動・共変動行列
	b <- (var(data)*(n-1))-w						# 群間変動・共変動行列
	within.cov <- w/(n-k)							# プールされた分散・共分散行列
	sd <- sqrt(diag(within.cov))
	within.r <- within.cov/outer(sd, sd)
	result <- geneig(b, w)							# 一般化固有値問題を解く
	eigen.values <- result$values						# 固有値
	nax <- sum(eigen.values > 1e-10)					# 解の個数
	eigen.values <- eigen.values[1:nax]					# 必要な個数の固有値を取り出す
	eigen.vectors <- result$vectors[, 1:nax,drop=FALSE]			# 固有ベクトルも取り出す
	LAMBDA <- rev(cumprod(1/(1+rev(eigen.values))))				# Wilks のΛ
	chi.sq <- ((p+k)/2-n+1)*log(LAMBDA)					# カイ二乗値
	l.L <- 1:nax-1
	df <- (p-l.L)*(k-l.L-1)							# 自由度
	p.wilks <- pchisq(chi.sq, df, lower.tail=FALSE)				# P 値
	canonical.corr.coef <- sqrt(eigen.values/(1+eigen.values))		# 正準相関係数
	temp <- diag(t(eigen.vectors) %*% w %*% eigen.vectors)
	temp <- 1/sqrt(temp/(n-k))
	coeff <- eigen.vectors*temp						# 正準判別係数
	const <- as.vector(-grand.means %*% coeff)				# 定数項
	coefficient <- rbind(coeff,const)					# 両方併せる
	std.coefficient <- coeff * sqrt(diag(within.cov))			# 標準化正準判別係数
	centroids <- t(t(coeff) %*% group.means+const)				# 重心
	can.score <- t(t(data %*% coeff)+const)					# 正準得点
	structure <- matrix(0, p, nax)						# 構造行列
	for (i in 1:nax) {
		ss <- split(data.frame(cs=can.score[,i], data), group)
		ss <- Reduce("+", lapply(ss, function(x) var(x)*(nrow(x)-1)))
		structure[,i] <- (ss[,1] / sqrt(ss[1,1] * diag(ss)))[-1]
	}
	d <- sapply(1:k, function(i) colSums((t(can.score)-centroids[i,])^2))	# 二乗距離
	p <- pchisq(d, p, lower.tail=FALSE)					# ケースがそれぞれの群に属するとしたとき，その判別値を取る確率
	p.Bayes <- t(t(exp(-d/2))*as.numeric(n.i))				# 各ケースが各群に所属するベイズ確率
	p.Bayes <- p.Bayes/rowSums(p.Bayes)					# P 値
	gname <- levels(group)
	classification <- factor(gname[max.col(p.Bayes)], levels=gname) 
	colnames(means) <- c("grand mean", gname)
	colnames(univariate) <- c("Wilks のラムダ", "F 値", "df1", "df2", "P 値")
	rownames(std.coefficient) <- rownames(structure) <- rownames(univariate) <- rownames(means)
	colnames(std.coefficient) <- colnames(structure) <- colnames(coefficient) <- 
	  colnames(centroids) <- colnames(can.score) <- paste("axis", 1:nax)
	rownames(coefficient) <- c(rownames(means), "constant")
	rownames(centroids) <- colnames(p.Bayes) <- colnames(p) <- gname
	rownames(can.score) <- rownames(p.Bayes) <- rownames(p) <- paste("case", 1:n)
	return(structure(list(means=means, univariate=univariate, betwee.ss=b, within.ss=w, within.cov=within.cov,
			within.r=within.r, eigen.values=eigen.values, LAMBDA=LAMBDA,
			chi.sq=chi.sq, df=df, p.wilks=p.wilks,
			canonical.corr.coef=canonical.corr.coef, std.coefficient=std.coefficient,
			structure=structure, coeff=coefficient, centroids=centroids, can.score=can.score,
			p.Bayes=p.Bayes, p=p, classification=classification, group=group, ngroup=k, nax=nax),
			class="candis"))
}
# print メソッド
print.candis <- function(	obj,						# candis が返すオブジェクト
				digits=5)					# 結果の表示桁数
{
	cat("\nWilks のラムダ\n\n")
	d <- data.frame(1:obj$nax, obj$LAMBDA, obj$chi.sq, obj$df, format.pval(obj$p.wilks))
	colnames(d) <- c("関数", "Wilks のラムダ", "カイ二乗値", "自由度", "P 値")
	print(d, digits=digits, row.names=FALSE)
	cat("\n判別係数\n\n")
	print(round(obj$coeff, digits=digits))
	cat("\n標準化判別係数\n\n")
	print(round(obj$std.coefficient, digits=digits))
	cat("\n構造行列\n\n")
	print(round(obj$structure, digits=digits))
	cat("\n判別結果\n\n")
	print(xtabs(~obj$group+obj$classification))
}
# summary メソッド
summary.candis <- function(	obj,						# candis が返すオブジェクト
				digits=5)					# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.candis <- function(obj,							# candis が返すオブジェクト
			pch=1:obj$ngroup,					# 3 群以上の判別時の記号
			col=1:obj$ngroup,					# 3 群以上の判別時の色
			xpos="topright", ypos=NULL,				# 3 群以上の判別時の legend 関数の x, y 引数
			which=c("boxplot", "barplot"),				# 2 群判別のときのグラフの種類
			nclass=20,						# barplot のおよその階級数
			...)							# plot, boxplot, barplot への引数
{
	score <- obj$can.score
	group <- obj$group
	if (ncol(score) >= 2) {
		int.group <- as.integer(group)
		plot(score, pch=pch[int.group], col=col[int.group], ...)
		legend(x=xpos, y=ypos, legend=levels(group), pch=pch, col=col)
	}
	else {
		which <- match.arg(which)
		if (which == "boxplot") {			# boxplot
			plot(score ~ group, xlab="群", ylab="判別値", ...)
		}
		else { 						# barplot
			tbl <- table(group, cut(score, breaks=pretty(score, n=nclass)))
			barplot(tbl, beside=TRUE, legend=TRUE, ...)
		}
	}
}
# Cutler - Ederer 法による生命表
# 原データから作成
ce.surv<- function(	time,					# 生存期間ベクトル
			event,					# 死亡なら 1，生存なら 0 の値をとるベクトル
			width)					# 生存率を計算する区間幅
{
	OK <- complete.cases(time, event)			# 欠損値を持つケースを除く
	time <- time[OK]
	event <- event[OK]
	ni <- length(time)					# 初期例数
	time <- floor(time/width)				# 生存期間を区間に分ける
	tbl <- table(factor(time, level=0:max(time)), event)	# 集計して二次データを作り，下請け関数を呼ぶ
	ce.surv2(ni, tbl[,2], tbl[,1], seq(0, width*(nrow(tbl)-1), width))
	
}
# 二次データから作成
ce.surv2 <- function(	ni,					# 初期例数
			d,					# 死亡数ベクトル
			u,					# 打ち切り数ベクトル
			interval)				# 区間数値ベクトル
{
	k <- length(d)			
	stopifnot(length(u) == k && length(interval) == k)	# 長さが同じであること
	n <- rep(ni, k)-cumsum(d+u)				# 各区間の開始時における例数
	n <- c(ni, n[-k])
	np <- n-u/2						# 区間中央の例数
	q <- d/np						# 死亡率
	p <- 1-q						# 生存率
	P <- cumprod(p)						# 累積生存率
	SE <- P*sqrt(cumsum(q/(n-u/2-d)))			# 累積生存率の標準誤差
	result <- data.frame(interval, n, d, u, np, q, p, P, SE)
	interval <- c(interval, max(interval)+interval[2]-interval[1])
	P <- c(1, P)
	plot(interval, P, type="o", xlim=c(0, max(interval)), ylim=c(0, 1))
	return(result)
}
# 検証的因子分析（確認的因子分析）
cfa <- function(r,						# 相関係数行列
		n,						# サンプルサイズ
		loc,						# 因子負荷量または設定情報
		lim=0.999,					# 推定パラメータの値域 (-lim, lim)
		init.val=0.5) {					# 推定パラメータの初期値
	tr <- function(x) sum(diag(x))				# トレース
	get.chi.sq <- function(par, scaler=TRUE)		# optim 用の関数
	{
		loadings[loadings != 0] <- par[1:p1]		# 因子負荷量行列
		fac.cor <- diag(nfac)				# 因子間相関行列
		fac.cor[lower.tri(fac.cor)] <- par[-(1:p1)]
		fac.cor <- fac.cor+t(fac.cor)
		diag(fac.cor) <- 1
		rho <- loadings%*%fac.cor%*%t(loadings)		# 母相関係数行列
		diag(rho) <- 1
		temp <- det(rho)				# 行列式が負の値になるときの対策
		retval <- tr(solve(rho)%*%r)+			# カイ二乗値の一部（定数部を除く）
		          (if (temp <= 0) 100 else log(temp))	# 100 は根拠のない仮の値
		if (scaler) {
			return(retval)				# optim への戻り値
		}
		else {
			return(list(retval=retval,		# 収束後の結果諸々
			            loadings=loadings, fac.cor=fac.cor, rho=rho))
		}
	}
	p <- nrow(r)						# 変数の個数
	if (is.matrix(loc) && nrow(loc) == p) {			# 因子負荷量として与えられるとき
		loadings <- loc					# 変数が因子に含まれるところを 1 とする
		nfac <- ncol(loadings)
	}
	else {
		nfac <- max(loc)				# 因子の個数
		loadings <- matrix(0, p, nfac)			# 因子負荷量行列
		loadings[cbind(1:p, loc)] <- 1			# 変数が因子に含まれるときに 1 とする
	}
	p1 <- sum(loadings != 0)				# 推定すべき因子負荷量の個数
	if (length(init.val) == 1) {				# 推定すべき因子負荷量の初期値
		init.val <- rep(init.val, p1)			# 全て同じにする
	}
	par <- c(init.val, rep(0, nfac*(nfac-1)/2))		# パラメータ初期値（因子負荷量行列と因子間相関行列）
	df <- p*(p-1)/2-length(par)				# 自由度
	ans <- optim(par, get.chi.sq, method="L-BFGS-B",	# 最適化（-1 ～ 1 の制約付きで）
	             lower=rep(-lim, length(par)), upper=rep(lim, length(par)))
	if (ans$convergence) {
		stop(paste("convergence =", ans$convergence,	# 何らかの原因で収束しなかったら，理由を明示し停止
	                   "\nmessage =", ans$message))
	}
	ans2 <- get.chi.sq(ans$par, scaler=FALSE)		# 最適パラメータのときの結果
	rho <- ans2$rho						# 母相関係数行列
	chisq <- (ans2$retval-log(det(r))-p)*(n-1)		# 定数部を含めたカイ二乗値
	P <- pchisq(chisq, df, lower.tail=FALSE)		# P 値
	z1 <- solve(rho)%*%(r-rho)
	z2 <- solve(rho)%*%r
	GFI <- 1-tr(z1%*%z1)/tr(z2%*%z2)			# GFI
	return(list(loadings=ans2$loadings, fac.cor=ans2$fac.cor,
	            chisq=chisq, P=P, df=df, GFI=GFI,
	            AGFI=1-p*(p+1)*(1-GFI)/2/df,		# AGFI
	            SRMR=sqrt(sum((rho-r)^2)/p/(p+1)),		# SRMR
	            RMSEA=sqrt((chisq-df)/df/(n-1)) ))	# RMSEA
}
# チョウ検定
chow <-function(dat1,					# 第 1 データセット
		dat2)					# 第 2 データセット
{
	ess <- function(dat)				# 残差平方和を返す関数
	{
		nc <- ncol(dat)
		ans <- lm(dat[,nc] ~ dat[,-nc])		# 右端の列が従属変数
		return(sum(ans$residuals^2))		# 残差平方和
	}

	method <- "Chow  検定"
	data.name <- paste(deparse(substitute(dat1)), "and", deparse(substitute(dat2)))
	dat1 <- subset(dat1, complete.cases(dat1))	# 欠損値を持つケースを除く
	dat2 <- subset(dat2, complete.cases(dat2))	# 欠損値を持つケースを除く
	ess12 <- ess(dat1)+ess(dat2)
	essc <- ess(rbind(dat1, dat2))
	df1 <- ncol(dat1)				# 第 1 自由度
	df2 <- nrow(dat1)+nrow(dat2)-2*df1		# 第 2 自由度
	f <- (essc-ess12)*df2/(df1*ess12)		# 検定統計量
	p <- pf(f, df1, df2, lower.tail=FALSE)		# P 値
	return(structure(list(statistic=c(F=f),
		parameter=c(df1=df1, df2=df2), p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# 母相関係数の信頼限界値を読みとるノモグラム
cl.r <- function(	conf,				# 信頼度（%で入力）
			n=c(5, 10, 30, 100, 500))	# 標本サイズ
{
	r2 <- seq(-1, 1, 0.2)				# 刻み幅 0.2
	r <- seq(-1, 1, 0.05)				# 刻み幅 0.05
	z <- atanh(r)					# Fisher の Z 変換値
	z0 <- atanh(0)					# 文字を書く位置を決めるため

	plot(r, r, type="n", xaxt="n", yaxt="n",	# 位置決めなど
	     xlab="r", ylab="rho",
	     main=paste(conf, "% confidence interval"))
	abline(h=r, v=r, col="pink")			# 粗い格子を描く
	abline(h=r2, v=r2, col="pink3")			# 細かい格子を描く

	sapply(n,					# 標本サイズごとに信頼曲線を描く
	       function(n)
	       {
		temp <- qnorm(0.5+conf/200)/sqrt(n-3)
		lines(r, tanh(z+temp), col="blue")
		lines(r, tanh(z-temp), col="blue")
		text(0, tanh(z0+temp), paste("n =", n), pos=2)
		text(0, tanh(z0-temp), paste("n =", n), pos=4)
	       })
	axis(1, r2)					# 横軸目盛り
	axis(2, r2)					# 縦軸目盛り
}
# データを，ほぼ同数のデータを含む 5 区分に区切り，塗り分け地図を描く
color.map1 <- function(	x,						# 長さ 47 の，統計データのベクトル
			color.no=8)					# 塗りつぶしに，何色系統を使うか（以下を参照）
{
	color.set <- matrix(c(						# 色の系統　color.no
		"gray100", "gray75", "gray50", "gray25", "gray0",	# 灰色1		1
		"#eeeeee", "#bbbbbb", "#999999", "#777777", "#555555",	# 灰色2		2
		"#ee0000", "#bb0000", "#990000", "#770000", "#550000",	# 赤色系	3
		"#00ee00", "#00bb00", "#009900", "#007700", "#005500",	# 緑色系	4
		"#0000ee", "#0000bb", "#000099", "#000077", "#000055",	# 青色系	5
		"#ee00ee", "#bb00bb", "#990099", "#770077", "#550055",	# 紫色系	6
		"#00eeee", "#00bbbb", "#009999", "#007777", "#005555",	# 水色系	7
		"#eeee00", "#bbbb00", "#999900", "#777700", "#555500"	# 黄色系	8
	), byrow=TRUE, ncol=5)

	if (!(color.no %in% 1:8)) {
		stop("color.no は，1～8 の整数でなければなりません")
	}
	div <- c(9, 19, 28, 38)
	xs <- sort(x)
	div2 <- (xs[div]+xs[div+1])/2
	map(1:47, color=color.set[color.no, findInterval(x, div2)+1])
}
# 統計データを 5 段階に区切って，塗り分け地図を描く
color.map2 <- function(	x,						# 長さ 47 の，統計データのベクトル
			t,						# データを 5 区分するための 4 つの値
			color.no=8)					# 塗りつぶしに，何色系統を使うか（以下を参照）
{
	color.set <- matrix(c(						# 色の系統　color.no
		"gray100", "gray75", "gray50", "gray25", "gray0",	# 灰色1		1
		"#eeeeee", "#bbbbbb", "#999999", "#777777", "#555555",	# 灰色2		2
		"#ee0000", "#bb0000", "#990000", "#770000", "#550000",	# 赤色系	3
		"#00ee00", "#00bb00", "#009900", "#007700", "#005500",	# 緑色系	4
		"#0000ee", "#0000bb", "#000099", "#000077", "#000055",	# 青色系	5
		"#ee00ee", "#bb00bb", "#990099", "#770077", "#550055",	# 紫色系	6
		"#00eeee", "#00bbbb", "#009999", "#007777", "#005555",	# 水色系	7
		"#eeee00", "#bbbb00", "#999900", "#777700", "#555500"	# 黄色系	8
	), byrow=TRUE, ncol=5)

	if (length(t) != 4) {
		stop("t は，長さ4のベクトルでなければなりません")
	}
	if (!(color.no %in% 1:8)) {
		stop("color.no は，1～8 の整数でなければなりません")
	}
	map(1:47, color=color.set[color.no, findInterval(x, t)+1])
}
# データにより塗り分け地図を描く
color.map3 <- function(	x,		# 長さ 47 の統計データベクトル
			t,		# データを区分する値
			color.set)	# 各区分を塗る色
{
	if (length(t)+1 != length(color.set)) {
		stop("t の長さは color.set の長さより 1 だけ小さくなければならない")
	}
	map(1:47, color=color.set[findInterval(x, t)+1])
}
# いくつかの都道府県を選択して描画，選択的に色づけも
color.map4 <- function(	prefs,			# 描画する都道府県番号のベクトル
			marks,			# 色づけをする都道府県のベクトル（prefs の部分集合）
			color,			# marks を塗る色
			others = "white")	# marks 以外の都道府県に塗る色
{
	map(prefs, color=ifelse(prefs %in% marks, color, others))
}
# 二次データに基づき，ピアソンの積率相関係数，スピアマンの順位相関係数，ケンドールの順位相関係数の無相関検定を行う
cor2.test <- function(	n,							# サンプルサイズ
			r,							# 相関係数
			conf.level=0.95,
			method = c("pearson", "kendall", "spearman"))		# 相関係数の種類
{
	data.name <- sprintf("n = %s, r = %s", n, r)
	method <- match.arg(method)						# 引数の補完
	if (method != "kendall") {						# ケンドールの順位相関係数以外の場合
		method <- paste(if (method == "pearson") "ピアソンの積率相関係数"
				else "スピアマンの順位相関係数", "の検定と推定", sep="")
		t <- abs(r)*sqrt((n-2)/(1-r^2))					# 検定統計量
		df <- n-2							# 自由度
		p <- pt(t, df, lower.tail=FALSE)*2				# P 値
		q <- qnorm(0.5-conf.level/2)
		conf.int <- tanh(atanh(r)+c(q, -q)/sqrt(n-3))
		attr(conf.int, "conf.level") <- conf.level
		return(structure(list(statistic=c(t=t), parameter=c(df=df), p.value=p,
			conf.int=conf.int, method=method, data.name=data.name),
			class="htest"))						# 結果をまとめて返す
	}
	else {									# ケンドールの順位相関係数の場合
		method <- "ケンドールの順位相関係数の検定"
		z <- abs(r)/sqrt((4*n+10)/(9*n*(n-1)))				# 検定統計量
		p <- pnorm(z, lower.tail=FALSE)*2				# P 値
		return(structure(list(statistic=c("Z-value"=z), p.value=p,
			method=method, data.name=data.name), class="htest"))	# 結果をまとめて返す
	}
}
# 相関比と決定係数を求める
correlation.ratio <- function(	x,						# 変数ベクトル
				group)						# 群を表す変数ベクトル
{
	ok <- complete.cases(x, group)						# 欠損値を持つケースを除く
	x <- x[ok]
	group <- factor(group[ok])
	n.i <- tapply(x, group, length)						# 各群のデータの個数
	n <- sum(n.i)								# 全データ数
	v.i <- tapply(x, group, var)						# 各群の不偏分散
	R.sq <- 1-sum((n.i-1)*v.i)/var(x)/(n-1)					# 決定係数
	c.r <- sqrt(R.sq)							# 相関比
	return(c("correlation ratio"=c.r, "coefficient of determination(R^2)"=R.sq))
}
# 共分散分析
covar.test <- function(	dat,						# データ行列
			cp1,						# 独立変数の列番号
			cp2,						# 従属変数の列番号
			cp3)						# 群変数の列番号
{
	dat <- subset(dat, complete.cases(dat[,c(cp1, cp2, cp3)]))	# 欠損値を持つケースを除く
	x <- dat[,cp1]							# 独立変数
	y <- dat[,cp2]							# 従属変数
	g <- dat[,cp3]							# 群変数

	nj <- table(g)							# 各群の例数
	n <- sum(nj)							# 全例数
	k <- length(nj)							# 群の個数
								# 独立変数について
	mx <- mean(x)							# 全体の平均値
	mxj <- tapply(x, g, mean)					# 各群の平均値
	sstx <- (n-1)*var(x)						# 全変動
	ssbx <- sum(nj*(mxj-mx)^2)					# 群間変動
	sswx <- sstx-ssbx						# 群内変動
								# 従属変数について
	my <- mean(y)							# 全体の平均値
	myj <- tapply(y, g, mean)					# 各群の平均値
	ssty <- (n-1)*var(y)						# 全変動
	ssby <- sum(nj*(myj-my)^2)					# 群間変動
	sswy <- ssty-ssby						# 群内変動

	spt <- (n-1)*cov(x, y)						# 共変動
	spb <- sum(nj*(mxj-mx)*(myj-my))				# 群間共変動
	spw <- spt-spb							# 群内共変動
	
	ss.wy <- sswy-spw^2/sswx					# 全群に共通な傾きと各群ごとの切片を持つ回帰直線からの変動
	ss.ty <- ssty-spt^2/sstx					# 全データに基づく回帰直線からの変動
	ss.by <- ss.ty-ss.wy						# 各群の回帰の差に起因する変動

	hensa.x <- x-mxj[g]						# 各群の平均値からの偏差
	hensa.y <- y-myj[g]						# 各群の平均値からの偏差
	xy <- hensa.x*hensa.y 
	xx <- hensa.x^2
	numerator <- tapply(xy, g, sum)
	denominator <- tapply(xx, g, sum)
	a <- numerator/denominator					# 各群のデータに基づく回帰直線の傾き
	b <- myj-a*mxj							# 各群のデータに基づく回帰直線の切片
	predict.y <- a[g]*x+b[g]					# 各群ごとのデータに基づく回帰直線による予測値
	ss.wyj <- sum((y-predict.y)^2)					# 各群ごとのデータに基づく回帰直線からの変動

	df.r <- k-1
	df.e <- n-2*k
	df.t <- n-k-1

	diff.reg <- ss.wy-ss.wyj					# 各群の回帰の差による推定誤差平方和
	ms.r <- diff.reg/df.r						# 各群の回帰の差による推定誤差分散
	ms.e <- ss.wyj/df.e						# 各群の推定誤差の和による推定誤差分散
	ms.w <- ss.wy/df.t						# 平均回帰に基づく推定誤差による推定誤差分散
	f <- ms.r/ms.e							# 検定統計量
	p <- pf(f, df.r, df.e, lower.tail=FALSE)			# P 値
	anova.table <- data.frame(					# 分散分析表としてまとめる
				"SS"=c(diff.reg, ss.wyj, ss.wy),
				"d.f."=c(df.r, df.e, df.t),
				"MS"=c(ms.r, ms.e, ms.w)
				 )
	rownames(anova.table) <- c("group x slope", "residual", "total")# 行の名前
	test.result <- c("F value"=f, "d.f.1"=df.r, "d.f.2"=df.e, "P value"=p)

	if (p <= 0.05) {						# 各群における回帰直線が等しくないとき，
		part2 <- anova.table2 <- test.result2 <- NULL		# それ以上の検定は行わない
	}
	else {								# 各群における回帰直線の傾きが等しくないとはいえないとき
		part2 <- "H0: 共変量で調整した平均値は同じである"	# 帰無仮説
		ms.by <- ss.by/df.r
		ms.wy <- ss.wy/df.t
		anova.table2 <- data.frame(				# 分散分析表としてまとめる
				"SS"=c(ss.by, ss.wy, ss.ty),
				"d.f."=c(df.r, df.t, n-2),
				"MS"=c(ms.by, ms.wy, ss.ty/(n-2))
				 )
		rownames(anova.table2) <- c("effect & group", "residual", "total")# 行の名前
		f2 <- ms.by/ms.wy					# 検定統計量
		p2 <- pf(f2, df.r, df.t, lower.tail=FALSE)			# P 値
		test.result2 <- c("F value"=f2, "d.f.1"=df.r, "d.f.2"=df.t, "P value"=p2)
	}
	return(list(	part1="H0: 各群の回帰直線の傾きは同じである",
			result1.1=anova.table,
			result1.2=test.result,
			part2=part2,
			result2.1=anova.table2,
			result2.2=test.result2))
}
#####
#
# クロス集計表を作成し，独立性の検定または代表値の差の検定を行う
#
#####

cross <- function(i,								# 表側に来る変数が入っているデータフレーム上の列番号または変数名ベクトル
		  j,								# 表側に来る変数が入っているデータフレーム上の列番号または変数名ベクトル
		  df,								# データフレーム
		  row=TRUE,							# 行ごとに 100% となるようにパーセントを取る
		  latex=TRUE,							# LaTeX 形式で度数分布表を出力する（デフォルトは LaTeX 形式）
                  captions=NULL,						# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
                  labels=NULL,							# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
		  test=c("none", "chisq", "fisher", "kruskal"),			# デフォルト none では検定を行わない。検定を行うときはその種類を指定する
		  output="",							# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
		  encoding=getOption("encoding"))				# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

# 下請け関数
	cross.sub <- function(ii, jj)						# ii, jj はスカラー。ii, jj で指定されたクロス集計を 1 つだけ行う
	{
		tbl <- table(df[,ii], df[,jj])					# 表の本体は table 関数で作る
		tbl <- tbl[rowSums(tbl) > 0,, drop=FALSE]			# 行和が 0 になる行を除く（factor 関数の使い方によってはこのような集計表ができる）
		tbl <- tbl[,colSums(tbl) > 0, drop=FALSE]			# 列和が 0 になる列を除く（同上）
		ans <- addmargins(tbl)						# 周辺和を付け加える
		nr <- nrow(ans)							# 集計表の行数
		nc <- ncol(ans)							# 集計表の列数
		colnames(ans)[nc] <- rownames(ans)[nr] <- "合計"		# 表頭，表側の該当箇所を「合計」とする
		pct <- ans*100 / if (row) ans[,nc] else rep(ans[nr,], each=nr)	# row の指示により，行 % か列 % のいずれかを取る
		if (latex) {							# LaTeX 形式で集計結果を出力する
			cat("\n\\begin{table}[htbp]\n", file=output)						# \begin{table}[htbp]
			if (is.null(captions)) {
				cat(sprintf("\\caption{%s ： %s}\n", colnames(df)[ii], colnames(df)[jj]), file=output)	# \caption{変数名 : 変数名}
			}
			else {
				cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
			}
			if (!is.null(labels)) {
				cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
			}
			cat("\\centering\n", file=output)							# \centering
			cat("\\begin{tabular}{l", rep("c", nc), "} \\hline\n", sep="", file=output)		# \begin{tabular}{cc…c} \hline
			cat(sprintf("& \\multicolumn{%i}{c}{%s}\\\\ \\cline{2-%i}\n", nc-1, colnames(df)[jj], nc), file=output)
														# 表頭の変数名
			cat(colnames(df)[ii], colnames(ans), sep=" & ", file=output)				# 表側の変数名 & 表頭1 & 表頭2 & … & 合計
			cat("\\\\ \\hline\n", file=output)							# \\ \hline
			for (i in 1:nr) {									# 各行について，
				cat(rownames(ans)[i], ans[i,], sep=" & ", file=output)				# 表側i & 観察数i1 & 観察数i2 & … & 合計 
				cat("\\\\\n", file=output)							# \\
				cat("\\%", sprintf("{\\small \\textit{%.1f}}", pct[i,]), sep=" & ", file=output)# % & パーセントi1 & パーセントi2 & … & パーセント
				cat("\\\\", file=output)							# \\
				if (i >= nr-1) {
					cat("\\hline\n", file=output)						# \hline \n
				}
				else {
					cat("\n", file=output)							# そのまま改行 \n
				}
			}
			cat("\\end{tabular}\n", file=output)							# \end{tabular}
		}
		else {								# tab で区切って出力する
			cat("\n表　", colnames(df)[ii], "：", colnames(df)[jj], sep="", file=output)		# 表　変数名：変数名
			cat("\n", colnames(df)[jj], sep="\t", file=output, fill=TRUE)				# 表頭の変数名
			cat(colnames(df)[ii], colnames(ans), sep="\t", file=output, fill=TRUE)			# 表側の変数名	表頭1	表頭2	…	合計
			for (i in 1:nr) {									# 各行について
				cat(rownames(ans)[i], ans[i,], sep="\t", file=output, fill=TRUE)		# 表側i	観察数i1	観察数i2	…	合計
				cat("%", sprintf("%.1f", pct[i,]), sep="\t", file=output, fill=TRUE)		# %	パーセントi1	パーセントi2	…	パーセント
			}
		}
		if (nr > 2 && nc > 2 && test != "none") {			# 2 行× 2 列以上の集計表については，検定オプションあり
			if (latex) {						# LaTeX 形式の出力なら表の後に追加
				cat("\\\\ \\noindent\n", file=output)
			}
			if (test == "chisq") {					# 独立性の検定 chisq を選んだ場合
				res <- chisq.test(tbl)				# chisq.test を使う
				cat(sprintf(if (latex) "$\\chi^2$値 = %.3f, 自由度 = %i, $P$値 = %.3f\n"
					    else "カイ二乗値 = %.3f, 自由度 = %i, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
			else if (test == "fisher") {				# 独立性の検定（Fisher の正確検定） fisher を選んだ場合
				cat(sprintf(if (latex) "$P_{Fisher}$値 = %.3f\n"
					    else "P 値（Fisher）= %.3f\n",
					    fisher.test(tbl)$p.value), file=output)
			}
			else if (test == "kruskal") {				# クラスカル・ウォリスの検定 kruskal を選んだ場合
				if (row) {					# 行ごとの % が 100% となるようにした row=TRUE の場合
					if (nc > 3 && (!is.ordered(df[,jj]) && !is.numeric(df[,jj]))) {
						warning(paste("「", colnames(df)[jj], "」は，順序尺度・間隔尺度・比尺度変数でなくてはなりません。", sep=""))
					}
					res <- kruskal.test(rep(col(tbl), tbl), rep(row(tbl), tbl))
				}
				else {						# 列ごとの % が 100% となるようにした row=FALSE の場合
					if (nr > 3 && (!is.ordered(df[,ii]) && !is.numeric(df[,ii]))) {
						warning(paste("「", colnames(df)[ii], "」は，順序尺度・間隔尺度・比尺度変数でなくてはなりません。", sep=""))
					}
					res <- kruskal.test(rep(row(tbl), tbl), rep(col(tbl), tbl))
				}
				cat(sprintf(if (latex) "$\\chi^2_{kw}$値 = %.3f, 自由度 = %i, $P$値 = %.3f\n"
					    else "カイ二乗値(kw) = %.3f, 自由度 = %i, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
		}
		if (latex) {							# LaTeX 形式で集計結果を出力する場合は，
			cat("\\end{table}\n", file=output)			# \end{table}
		}
	}

	getNum <- function(str, df) {					# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

# cross 関数の本体
	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

	test <- match.arg(test)							# test 引数から，完全な検定手法名を得る
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}
	index <- 0
	for (ii in i) {								# i はベクトルまたはスカラー
		for (jj in j) {							# j はベクトルまたはスカラー
			if (ii != jj) {						# i, j の全ての組み合わせについて（ii と jj が違うときのみ），
				index <- index+1
				cross.sub(ii, jj)				# クロス集計のための下請け関数 cross.sub を呼ぶ
			}
		}
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
xtable.ftable <- function(	obj,				# crosstabsn が返す ftable オブジェクト
			caption="caption",			# キャプション
			label="label",			# ラベル
			percentage=c("row", "col", "none"),	# % を付ける方向
			same.line=TRUE, 			# % を度数と同じ行に付けるときは TRUE にする
			percentage.font=c("small", "footnotesize", "tiny", "normalsize"), # LaTeX でのフォントサイズの指定 tiny, footnotesize など
			position=c("c", "r", "l"), 		# フィールド内での配置 "c", "r", "l" のいずれか
			rev=-1.5,				# 行間を詰めるための，逆改行の大きさをミリ単位で指定（逆改行しない場合には 0 を指定する）
			booktabs=FALSE)			# TRUE なら \hline の代わりに \toprule, \midrule, \bottomrule を使う
# ftable 関数が返した ftable クラスのオブジェクトを入力し，LaTeX ソースを出力する
# formula で指定する。~ の左辺には1個のみ指定できる
# Sweave から使用するのが便利
# 使用例
# Titanic
# x <- ftable(Survived ~ ., data = Titanic)
# a <- ftable(Survived ~ Sex + Class + Age, data = x)
# xtable(a)
# xtable(ftable(Survived ~ Sex + Class, data = x))

{	row.vars <- attr(obj, "row.vars")
	n.row.vars <- length(row.vars)
	names.row.vars <- names(row.vars)
	m.row.vars <- sapply(row.vars, length)
	col.vars <- attr(obj, "col.vars")
	n.col.vars <- length(col.vars)
	names.col.vars <- names(col.vars)
	m.col.vars <- sapply(col.vars, length)
	if (n.col.vars != 1) {
		stop("col.vars が 1 変数の ftable オブジェクトしか扱えません")
	}
	nrow <- nrow(obj)
	side <- matrix("", nrow, n.row.vars)
	n.block <- nrow/m.row.vars[n.row.vars]
	side[, n.row.vars] <- unlist(rep(row.vars[n.row.vars], n.block))
	for (i in seq_len(n.row.vars-1)) {
		every <- prod(m.row.vars[(i+1):n.row.vars])
		side[(0:(nrow-1))%%every==0, i] <- unlist(row.vars[i])
	}

	percentage <- match.arg(percentage)
	if (percentage == "none") {
		same.line <- FALSE
	}
	percentage.font <- match.arg(percentage.font)
	position <- match.arg(position)

	if (booktabs) {
		toprule <- "\\toprule"
		midrule <- "\\midrule"
	}
	else {
		toprule <- midrule <- "\\hline"
	}

	col.vars <- c(unlist(col.vars[[1]]), "合計")
	fac <- same.line+1
	if (same.line) {
		pos <- c(rep(position, n.row.vars), rep(paste(position, "@{}", position), m.col.vars+1))
		header <- paste(paste(names.row.vars, collapse=" & "), paste("&", paste(col.vars, "\\%", sep=" & ", collapse=" & ")))
		fmt <- sprintf("%%d & {\\%s \\textit{%%6.1f}}", percentage.font)
	}
	else {
		pos <- rep(position, m.col.vars+1+n.row.vars)
		header <- paste(paste(names.row.vars, collapse=" & "), paste(col.vars, collapse=" & "), sep=" & ")
		fmt <- sprintf("{\\%s \\textit{%%5.1f}}", percentage.font)
	}
	cat("\\begin{table}[htbp]\n",
	    "\\caption{", caption, "}\n",
	    "\\label{", label, "}\n",
	    "\\centering\n",
	    "\\begin{tabular}{", pos, "} ", toprule, " \n", sep="")
	cat(paste(rep("&", n.row.vars), collapse=" "))
	cat(sprintf(" \\multicolumn{%i}{c}{%s}\\\\ \\cline{%i-%i}\n",
	    fac*m.col.vars[1], names.col.vars[1], n.row.vars+1, fac*m.col.vars[1]+n.row.vars))
	cat(header, " \\\\ ", midrule, "\n", sep="")

	for (k in 1:n.block) {
		end <- k*m.row.vars[n.row.vars]
		begin <- end-m.row.vars[n.row.vars]+1
		block <- addmargins(obj[begin:end, ])
		side.block <- rbind(side[begin:end, , drop=FALSE ], c(rep("", n.row.vars-1), "合計"))
		if (percentage == "row") {
			pct <- block/block[, m.col.vars+1]*100
		}
		else {
			pct <- t(t(block)/block[m.row.vars[n.row.vars]+1,]*100)
		}
		n <- m.row.vars[n.row.vars]+1
		for (i in 1:n) {
			cat(sprintf("%s &", side.block[i,]))
			if (same.line) {
				cat(gsub("NaN", "---", paste(apply(cbind(block[i,], pct[i,]), 1, function(y) sprintf(fmt, y[1], y[2])), collapse=" & ")))
			}
			else {
				cat(paste(block[i,], collapse=" & "), "\\\\ \n")
				if (percentage != "none") {
					cat(rep(" &", n.row.vars-1))
					cat("\\%", gsub("NaN", "---", sprintf(fmt, pct[i, ])), sep=" & ")
				}
			}
			if (percentage != "none") {
				cat(" \\\\")
			}
			if (i < n-1) {
				cat(sprintf("[%smm]\n", rev))
			}
			else if (i == n) {
				if (end < nrow) {
					cat(sprintf("\\cline{%i-%i}\n", sum(side[end+1,] == "")+1, fac*(m.col.vars[1]+1)+n.row.vars))
				}
				else {
					cat(sprintf("%s\n", toprule))
				}
			}
			else {
				cat(sprintf("\\cline{%i-%i}\n", n.row.vars, fac*(m.col.vars[1]+1)+n.row.vars))
			}
		}
	}
	cat("\\end{tabular}\n",
	    "\\end{table}\n", sep="")
}
# カルダーノ法により，3次方程式の解を求める
cubic <- function(a, b, c, d)
{
	cubic.root <- function(x) return(sign(x)*abs(x)^(1/3))			# a x^3 + b x^2 + c x +d = 0 の係数
	res <- NULL
	res$coefficients <- c(a, b, c, d)
	if (a == 0) {
		return("3次の項の係数がゼロです")
	}
	b <- b/(3*a)
	c <- c/a
	d <- d/a
	p <- b^2-c/3
	q <- (b*c-2*b^3-d)/2
	a <- q^2-p^3
	if (a == 0) {
		q <- cubic.root(q)
		res$ans <- c(2*q-b, -q-b)
	}
	else if (a > 0) {
		a3 <- cubic.root(q+sign(q)*sqrt(a))
		b3 <- p/a3
		x <- complex(real=-(a3+b3)/2-b, imaginary=abs(a3-b3)*sqrt(3)/2)
		res$ans <- c(a3+b3-b, x, Conj(x))
	}
	else {
		a <- 2*sqrt(p)
		t <- acos(q/(p*a/2))
		res$ans <- c(a*cos(t/3)-b, a*cos((t+2*pi)/3)-b, a*cos((t+4*pi)/3)-b)
	}
	class(res) <- "cubic"
	return(res)
}
# cubic クラスのオブジェクトを表示する
print.cubic <- function(obj)
{
	put0 <- function(x) return(paste(ifelse(x >= 0, "", "-"), ifelse(abs(x) == 1, "", abs(x)), sep=""))
	put1 <- function(x) return(paste(ifelse(x >= 0, "+", "-"), ifelse(abs(x) == 1, "", abs(x)), sep=""))
	put2 <- function(x) return(paste(ifelse(x >= 0, "+", "-"), abs(x), sep=""))
	printf("%sx^3%sx^2%sx%s=0\n", put0(obj$coefficients[1]), put1(obj$coefficients[2]),
	                              put1(obj$coefficients[3]), put2(obj$coefficients[4]))
	sapply(obj$ans, print)
}
# 行名（行番号）なしでデータフレーム（行列）を表示する
data.list <- function(d)		# データフレームまたは行列
{
	cat(paste(t(d), c(rep("\t", ncol(d)-1), "\n")), sep="")
}
# 行名（行番号）なしでデータフレーム（行列）を表示する
data.list2 <- function(d)		# データフレームまたは行列
{
	invisible(apply(d, 1, function(i) cat(paste(i, c(rep("\t", length(i)-1), "\n")))))
}
# 行名（行番号）なしでデータフレーム（行列）を表示する
data.list3 <- function(	d,		# データフレームまたは行列
			num = FALSE)	# 行番号をつけるかどうか
{
	nv <- ncol(d)
	n <- nrow(d)
	for (i in 1:n) {
		if (num) cat(i, "\t")
		for (j in 1:(nv-1)) {
			cat(d[i, j], "\t")
		}
		cat(d[i, j+1], "\n")
	}
}
# ユリウス日から何曜日かを求める
dw <- function(j)
{
	c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")[(j+1)%%7+1]
}
# 年月日からユリウス日を求める
# 以下の式で，整数演算（%/%）がかけ算（*）やわり算（/）より優先順位が高いことに注意
J.day <- function(iy, jm, kd)
{
	tmp <- -(jm < 3)
	kd-32075+(1461*(iy+4800+tmp))%/%4+(367*(jm-2-tmp*12))%/%12-(3*((iy+4900+tmp)%/%100))%/%4
}
# ユリウス日を年月日に直す
date2 <- function(jul)
{
	l <- jul+68569
	n <- (4*l)%/%146097
	l <- l-(146097*n+3)%/%4
	iy <- (4000*(l+1))%/%1461001
	l <- l-(1461*iy)%/%4+31
	jm <- (80*l)%/%2447
	kd <- l-(2447*jm)%/%80
	l <- jm%/%11
	jm <- jm+2-12*l
	iy <- 100*(n-49)+iy+l
	sprintf("%4g/%02g/%02g", iy, jm, kd)
}
# 同じサンプルからの二つの相関係数に差があるか検定する
diff.r <- function(	n,				# サンプルサイズ
			rxy, rvy,			# 差を計算する相関係数1
			rxv)				# 先の 2 つの相関係数に関するもう一つの相関係数
{
	method <- "同じサンプルからの二つの相関係数に差があるか検定する"
	data.name <- sprintf("n = %s, r_xy = %s, r_vy = %s, r_xv = %s", n, rxy, rvy, rxv)
	detR <- (1-rxy^2-rvy^2-rxv^2)+2*rxy*rxv*rvy	# 行列式
	t0 <- (abs(rxy-rvy)*sqrt((n-1)*(1+rxv)))/	# 検定統計量
		sqrt(2*detR*(n-1)/(n-3)+(rxy+rvy)^2*(1-rxv)^3/4)
	df <- n-3					# 自由度
	p <- pt(abs(t0), df, lower.tail=FALSE)*2	# P 値
	return(structure(list(statistic=c(t=t0),	# 結果をまとめて返す
		parameter=c(df=df), p.value=p, method=method,
		data.name=data.name), class="htest"))
}
# 対応のあるデータの二つの相関係数の相等性の検定
# http://koko15.hus.osaka-u.ac.jp/%7Ekano/lecture/faq/q1.html#excel
diff.r2 <- function(x)
{
	method <- "対応のあるデータの二つの相関係数の相等性の検定"
	data.name <- deparse(substitute(x))
	x <- as.matrix(x)
	n <- nrow(x)
	r <- cor(x)
	z12 <- atanh(r[1,2])
	z34 <- atanh(r[3,4])
	a1 <- r[1,3]*r[2,4]+r[1,4]*r[2,3]
	a2 <- -r[3,4]*(r[1,3]*r[2,3]+r[1,4]*r[2,4])
	a3 <- -r[1,2]*(r[1,3]*r[1,4]+r[2,3]*r[2,4])
	a4 <- r[1,2]*r[3,4]*(r[1,3]^2+r[1,4]^2+r[2,3]^2+r[2,4]^2)/2
	d <- (1-r[1,2]^2)*(1-r[3,4]^2)
	chi.sq <- (n-3)*(z12-z34)^2/(2-2*(a1+a2+a3+a4)/d)
	p <- pchisq(chi.sq, 1, lower.tail=FALSE)
	return(structure(list(statistic=c("chi sq."=chi.sq), parameter=c(df=1),
		p.value=p, method=method, data.name=data.name), class="htest"))
}
# 線形判別関数
disc <- function(	data,						# データ行列
			group,						# 群変数
			func.name=c("solve", "ginv"))			# 逆行列を計算する関数名
{
	inverse <- if (match.arg(func.name) == "solve") solve else { library(MASS); ginv}
	data <- as.data.frame(data)
	if (is.null(colnames(data))) {
		colnames(data) <- paste("Var", 1:p, sep="")
	}
	vname <- colnames(data)
	group <- as.factor(as.matrix(group))				# factor にする
	gname <- levels(group)
	OK <- complete.cases(data, group)				# 欠損値を持つケースを除く
	data <- as.matrix(data[OK,])
	group <- group[OK]
	p <- ncol(data)							# 説明変数の数
	ncase <- nrow(data)						# サンプルサイズ
	num <- table(group)						# 各群のサンプルサイズ
	ng <- length(num)						# 群の数
	g.name <- names(num)						# 群の名前
	means <- t(matrix(unlist(by(data, group, colMeans)), p))	
	g.mean <- colMeans(data)
	t <- var(data)*(ncase-1)					# 分散共分散行列
	vars <- by(data, group, function(x) var(x)*(nrow(x)-1))		# 各群の変動・共変動行列
	w <- matrix(rowSums(matrix(mapply("+", vars), ncol=ng)), p)	# 群内変動・共変動行列
	g.sd <- apply(data, 2, sd)					# 各変数の標準偏差
	det.w <- det(w)							# 行列式
	det.t <- det(t)							# 行列式
	wl <- det.w/det.t
	inv.w <- inverse(w)						# 逆行列
	a <- -2*(ncase-ng)*inv.w%*%t(means)
	a0 <- rowSums(means%*%inv.w*means)*(ncase-ng)
	c.function <- rbind(a, a0)
	temp <- matrix(0, nr=ng, nc=ng)
	temp1 <- row(temp)
	temp1 <- temp1[upper.tri(temp1)]
	temp2 <- col(temp)
	temp2 <- temp2[upper.tri(temp2)]
	d.function <- sapply(1:length(temp1), function(i) (c.function[,temp2[i]]-c.function[,temp1[i]])/2)
	F <- diag(inverse(t)/inv.w)
	idf1 <- ng-1							# 第一自由度
	idf2 <- ncase-idf1-p						# 第二自由度
	F <- idf2/idf1*(1-F)/F						# F 値
	P <- pf(F, idf1, idf2, lower.tail=FALSE)			# P 値
	c1 <- ifelse(p^2+idf1^2 != 5, sqrt((p^2*idf1^2-4)/(p^2+idf1^2-5)), 1)
	c2 <- wl^(1/c1)
	df1 <- p*idf1
	df2 <- (ncase-1-(p+ng)/2)*c1+1-0.5*p*idf1
	F.wl <- df2*(1-c2)/(df1*c2)
	P.wl <- pf(F.wl, df1, df2, lower.tail=FALSE)
	t.data <- t(data)
	D2 <- (ncase-ng)*matrix(sapply(1:ng, function(i) {temp <- t.data-means[i,]; sapply(1:ncase, function(j) temp[,j]%*%inv.w%*%temp[,j])}), nr=ncase)
	P2 <- pchisq(D2, p, lower.tail=FALSE)
	prediction <- as.factor(g.name[apply(P2, 1, order)[ng,]])
	correct <- ifelse(prediction == group, TRUE, FALSE)
	correct.table <- table(group, prediction)
	correct.rate <- sum(diag(correct.table))/ncase*100
	factor1 <- levels(group)
	factor2 <- levels(prediction)
	if (ng ==2) {							# 2 群の場合には，判別値を計算
		discriminant.value <- data%*%d.function[1:p]+d.function[p+1]
	}
	else {
		discriminant.value <- NULL
	}
	colnames(c.function) <- paste("g", 1:ng, sep="")
	colnames(D2) <- paste("D to", gname)
	colnames(P2) <- paste("P to", gname)
	colnames(d.function) <- paste(gname[temp1], gname[temp2], sep=":")
	rownames(c.function) <- rownames(d.function) <- c(vname, "constant")
	colnames(c.function) <- gname
	return(structure(list(d.function=d.function, c.function=c.function, partial.F=F, partial.F.P=P, df1=idf1, df2=idf2, wilks.lambda=wl, wilks.lambda.F=F.wl, wilks.lambda.P=P.wl, wilks.lambda.df1=df1, wilks.lambda.df2=df2, distance=D2, P.value=P2, prediction=prediction, correct=correct, correct.table=correct.table, correct.rate=correct.rate, discriminant.value=discriminant.value, group=group, factor1=factor1, factor2=factor2), class="disc"))
}
# print メソッド
print.disc <- function(	obj,						# disc 関数が返すオブジェクト
			digits=5)					# 結果の表示桁数
{
	cat("\n判別関数\n\n")
	result <- cbind(obj$d.function, "Partial F"=c(obj$partial.F, NA), "p-value"=c(obj$partial.F.P, NA))
	print.default(round(result, digits=digits), na.print="")
	cat("\n分類関数\n\n")
	print(round(obj$c.function, digits=digits))
	cat("\n判別結果\n\n")
	print(obj$correct.table)
	cat(sprintf("\n正判別率 = %.1f %%\n", obj$correct.rate))
}
# summary メソッド							# すべての結果を表示する
summary.disc <- function(	obj,					# disc が返すオブジェクト
				digits=5)				# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.disc <- function(	obj, 						# disc 関数が返すオブジェクトの
			which=c("boxplot", "barplot", "scatterplot"),	# 箱髭図か棒グラフか散布図かの選択
			nclass=20,					# 棒グラフの場合のおよその階級数
			pch=1:ncol(obj$distance),				# scatterplot を描く記号
			col=1:ncol(obj$distance),				# scatterplot の記号の色
			xpos="topright", ypos=NULL,			# scatterplot の凡例の位置
			...)						# boxplot, barplot, scatterplt に引き渡す引数
{
	if (!is.null(obj$discriminant.value)) {
		which <- match.arg(which)
		if (which == "boxplot") {				# boxplot
			plot(obj$discriminant.value ~ obj$group, xlab="群", ylab="判別値", ...)
		}
		else if (which == "barplot") { 				# barplot
			tbl <- table(obj$group, cut(obj$discriminant.value,
					breaks=pretty(obj$discriminant.value, n=nclass)))
			barplot(tbl, beside=TRUE, legend=TRUE, xlab="判別値", ...)
		}
		else {							# scatterplot 各群の重心までの二乗距離
			group <- obj$group
			group.levels <- levels(group)
			distance1 <- obj$distance[,1]
			distance2 <- obj$distance[,2]
			max1 <- max(distance1)
			max2 <- max(distance2)
			max0 <- max(max1, max2)
			plot(distance1, distance2, col=col[as.integer(group)], pch=pch[as.integer(group)],
				xlim=c(0, max0), xlab=paste(group.levels[1], "の重心への二乗距離"),
				ylim=c(0, max0), ylab=paste(group.levels[2], "の重心への二乗距離"), asp=1, ...)
			abline(0, 1, lty=2)
			text(max1, max2/2, paste(group.levels[2], "に判別される領域"), pos=2)
			text(0, max2+strheight("H")*1.5, paste(group.levels[1], "に判別される領域"), pos=4)
			legend(x=xpos, y=ypos, legend=group.levels, col=col, pch=pch)
		}
	}
	else {
		warning("3群以上の場合にはグラフ表示は用意されていません")
	}
}
# 整数の約数を見つける
divisor <- function(n)					# 整数
{
	if (n >= 2^53) return(NA)			# 限界を超えると NA を返す
	else if (n %% 2 == 0) return(2)
	else if(n %% 3 == 0) return(3)
	maxitr <- as.integer(sqrt(n))
	i <- 1
	repeat {
		i <- i+4
		if (i > maxitr) return(n)
		else if (n %% i == 0) return(i)
		i <- i+2
		if (i > maxitr) return(n)
		else if (n %% i == 0) return(i)
	}
}
# 素数の判定
is.prime <- function(n)
{
	return(ifelse(n >= 2^53, NA, n == divisor(n)))	# 限界を超えると NA を返す
}
# 素因子分解
factorization <- function(n, simple=TRUE)
{
	if (n >= 2^53) return(NA)			# 限界を超えると NA を返す
	result <- NULL
	repeat {
		div <- divisor(n)
		result <- append(result, div)
		n <- n/div
		if (n == 1) break
	}
	if (simple) {
		return(result)
	}
	else {
		result <- table(result)
		result <- paste(as.integer(names(result)), "^", result, collapse=" * ")
		result <- gsub("\\^ 1 ", "", result) 
		result <- sub(" \\^ 1$", "", result) 
		return(result)
	} 
}
# 度数分布表を作成し，ヒストグラムを描く
dosuu.bunpu <- function(x,				# データベクトル
			w,				# 階級幅
			percent=FALSE)			# TRUE にすると，縦軸が % 目盛りになる
{
	x <- x[!is.na(x)]				# 欠損値を持つケースを除く
	y <- floor(x/w)					# 階級に分ける
	mn <- min(y)					# 最小の階級
	mx <- max(y)					# 最大の階級
	y <- y-mn+1					# 最小値が1になるように変換
	freq <- table(factor(y, levels=1:(mx-mn+1)))	# 度数分布表（度数が 0 になる階級も確保）
	names(freq) <- mn:mx*w				# 階級名
	pcnt <- freq/sum(freq)*100			# パーセント
	cum.pcnt <- cumsum(pcnt)			# 累積パーセント
	h <- if (percent) freq else pcnt		# 縦軸の選択
	barplot(h, axis.lty="solid", space=c(0, 0))	# ヒストグラムとして描く
	return(cbind(freq, pcnt, cum.pcnt))
}
# 群別のデータプロット
dot.plot <- function(	x,						# 群変数ベクトル
			y,						# データベクトル
			accu=0,						# データを階級化するための値
			stp=0,						# 水平方向に記号を並べるときのずらす量
			log.flag=FALSE,					# 縦軸を対数目盛りにするとき TRUE
			simple=FALSE,					# 対数目盛りのとき，目盛数値を 10 のべき乗に限るなら TRUE
			symmetrical=TRUE,				# 記号を左右対称にするなら TRUE
			...)						# plot 関数に引き渡すその他の引数
{
	OK <- complete.cases(x, y)					# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	x.name <- unique(x)						# 群を表す変数の取る値（factor でありうる）
	if (is.factor(x)) {						# factor なら，
		x <- as.integer(x)					# 整数値に戻す
	}
	if (log.flag == TRUE) {						# 対数目盛りで描くなら，
		y0 <- y							# 値のバックアップをとってから，
		y <- log10(y)						# 常用対数をとる
	}
	if (accu == 0) {						# accu のデフォルト値を計算
		accu <- diff(range(y))/100				# 最大値と最小値の差の百分の一
	}
	if(stp == 0) {							# spt のデフォルト値を計算
		stp <- diff(range(x))/100				# 最大値と最小値の差の百分の一
	}
	y <- round(y/accu)*accu						# y を丸める
	x1 <- unique(x)							# 群を表す変数の種類（数値に変換したもの）
	for (i in seq(along=x1)) {					# 全ての群について，
		freq <- table(y[x==x1[i]])				# ある群のデータについて度数分布を求める
		for (j in seq(along=freq)) {				# 度数分布の各階級について
			if (freq[j] >= 2) {				# 複数個のデータがあるならば，
				offset <- ifelse(symmetrical, (freq[j]-1)/2*stp, 0)	# 対称に描くかどうかで描き始めが違う
				for (k in seq(along=y)) {
					if (abs(y[k]-as.numeric(names(freq)[j])) < 1e-10 && abs(x[k]-x1[i]) < 1e-10) {
						freq[j] <- freq[j]-1
						x[k] <- x[k]-offset+freq[j]*stp
					}
				}
			}
		}
	}
	if (log.flag) {							# 対数目盛りで描くなら，
		plot(x, y, type="n", xaxt="n", yaxt="n", ...)
		options(warn=-1)
		points(x, y, ...)
		options(warn=0)
		y0 <- floor(log10(y0))
		log.min <- min(y0)
		y2 <- 1:10*10^log.min
		n <- max(y0)-log.min
		y1 <- rep(y2, n+1)*10^rep(0:n, each=10)
		if (simple) {
			y2 <- y1[abs(log10(y1)-round(log10(y1))) < 1e-6]
			axis(2, at=log10(y1), labels=FALSE)
			axis(2, at=log10(y2), labels=y2)
		}
		else {
			axis(2, at=log10(y1), labels=y1)
		}
	}
	else {
		plot(x, y, xaxt="n", ...)
	}
	axis(1, at=x1, labels=as.character(x.name))
	print(paste("accu =", accu, "   stp = ", stp), quote=FALSE)
}
# 境界線データに基づいて白地図を描く
draw.map <- function(fn)							# 境界線データのあるファイル名
{
	data <-read.table(fn, header=FALSE) 					# x, y 座標が組みになっている
	data[data[,1]==0 & data[,2]==0,] <- NA					# x, y 座標が共に 0 であるのは，一連の境界線の終わりを意味する
	plot(data, type = "l", axes=FALSE, bty="n", xlab="", ylab="", asp=1)	# (NA, NA) は，結線されない
}
# 二次元クロス集計表・度数表を双対尺度法で分析する
dual <- function(tbl)						# 二次元クロス表（度数表）を行列として与える
{
	tbl <- data.matrix(tbl)					# データフレームも行列にする
	nr <- nrow(tbl)						# 行数
	nc <- ncol(tbl)						# 列数
	if (is.null(rownames(tbl))) {				# 行ラベルの補完
		rownames(tbl) <- paste("Row", 1:nr, sep="-")
	}
	if (is.null(colnames(tbl))) {				# 列ラベルの補完
		colnames(tbl) <- paste("Col", 1:nc, sep="-")
	}
	ft <- sum(tbl)						# サンプルサイズ
	fr <- rowSums(tbl)					# 行和
	fc <- colSums(tbl)					# 列和
	temp <- sqrt(outer(fc, fc))
	w <- t(tbl/fr)%*%tbl/temp-temp/ft
	res <- eigen(w)						# 固有値・固有ベクトルを求める
	ne <- length(res$values[res$values > 1e-5])		# 有効な解の個数
	vec <- res$vectors[,1:ne, drop=FALSE]			# 固有ベクトルの切り詰め
	val <- res$values[1:ne]					# 固有値の切り詰め
	col.score <- vec*sqrt(ft/fc)				# 列スコア
	row.score <- tbl%*%col.score/outer(fr, sqrt(val))	# 行スコア
	col.score2 <- t(t(col.score)*sqrt(val))			# 相関比で重み付けした列スコア
	row.score2 <- t(t(row.score)*sqrt(val))			# 相関比で重み付けした列スコア
	cont <- val/sum(diag(w))*100				# 寄与率
	chi.sq <- -(ft-1-(nr+nc-1)/2)*log(1-val)		# カイ二乗値
	df <- nr+nc-1-2*(1:ne)					# 自由度
	P <- pchisq(chi.sq, df, lower.tail=FALSE)		# P 値
	result <- matrix(c(val, sqrt(val), cont, cumsum(cont), chi.sq, df, P), byrow=TRUE, ncol=ne)
	rownames(result) <- c("eta square", "correlation", "contribution", "cumulative contribution", "Chi square value", "degrees of freedom", "P value")
	colnames(result) <- colnames(row.score) <- colnames(col.score) <- paste("Axis", 1:ne, sep="-")
	rownames(row.score) <- rownames(tbl)
	rownames(col.score) <- colnames(tbl)
	dimnames(row.score2) <- dimnames(row.score)
	dimnames(col.score2) <- dimnames(col.score)
	result <- list(result=result, row.score=row.score, column.score=col.score, row.score.weighted=row.score2, column.score.weighted=col.score2)
	class(result) <- c("dual")
	invisible(result)
}
# ダネットの方法による多重比較
dunnett <- function(	data,					# データベクトル
			group)					# 群変数ベクトル
{
	get.rho <- function(ni)					# ρを計算する
	{
		k <- length(ni)
		rho <- outer(ni, ni, function(x, y) { sqrt(x/(x+ni[1])*y/(y+ni[1])) })
		diag(rho) <- 0
		sum(rho[-1, -1])/(k-2)/(k-1)
	}

	pdunnett <- function(x, a, df, r)			# P 値を計算する
	{
		corr <- diag(a-1)
		corr[lower.tri(corr)] <- r
		1-pmvt(lower=-x, upper=x, delta=numeric(a-1), df=df, corr=corr, abseps=0.0001)
	}

	OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
	data <- data[OK]
	group <- factor(group[OK])				# 群変数は factor に変換
	ni <- table(group)					# 各群のデータ数
	a <- length(ni)						# 群の数
	n <- length(data)					# 全体のデータ数
	mi <- tapply(data, group, mean)				# 各群の平均値
	vi <- tapply(data, group, var)				# 各群の不偏分散
	Vw <- sum(vi*(ni-1))/(n-a)				# 群内分散
	rho <- get.rho(ni)					# ρ
	t <- (abs(mi-mi[1])/sqrt(Vw*(1/ni+1/ni[1])))[2:a]	# 対照群と各群の比較における t 値
	p <- sapply(t, pdunnett, a, n-a, rho)			# P 値
	result <- cbind(t, p)
	rownames(result) <- paste(1, 2:a, sep=":")
	return(result)
}
# ED50 や LD50 を求める
ed50 <- function(	x,					# 各群の用量ベクトル
			n,					# 各群の例数ベクトル
			r,					# 各群の反応数ベクトル
			transform=c("none", "ln", "log"))	# 対数変換をするかどうか
{
	# 第一近似
	estimation1 <- function(x, n, r)
	{
		select <- r > 0 & r < n				# 条件に合うものだけを取り出す
		x <- x[select]
		n <- n[select]
		r <- r[select]
		p9 <- r/n					# 標本比率
		y <- qnorm(p9)+5				# プロビット変換
		w <- 1
		wx2 <- sum(w*x^2)
		wy <- sum(w*y)
		wx <- sum(w*x)
		wyx <- sum(w*y*x)
		ww <- length(x)
		temp <- ww*wx2-wx^2
		alpha <- (wx2*wy-wx*wyx)/temp			# 切片
		beta <- (ww*wyx-wx*wy)/temp			# 傾き
		ollf <- sum(r*log(p9)+(n-r)*log(1-p9))		# 観察値の対数尤度
		list(ollf=ollf, alpha=alpha, beta=beta)
	}

	# 第二近似（第一近似から出発して収束計算）
	estimation2 <- function(x, n, r, alpha, beta, epsilon=1e-5)
	{
		sqrt.pi2 <- 1/sqrt(2*pi)
		for (ii in 1:100) {
			y <- alpha+beta*x
			x9 <- y-5
			p9 <- pnorm(x9)
			z <- exp(-x9^2/2)*sqrt.pi2
			y <- y+ (r/n-p9)/z
			w <- n*z^2/((1-p9)*p9)
			wx2 <- sum(w*x^2)
			wy <- sum(w*y)
			wx <- sum(w*x)
			wyx <- sum(w*y*x)
			ww <- sum(w)
			ssnwx <- wx2-wx^2/ww
			beta1 <- (wyx-wx*wy/ww)/ssnwx
			xbar <- wx/ww
			alpha1 <- (wy/ww)-beta1*xbar
			if (abs(alpha-alpha1)/alpha < epsilon && abs(beta-beta1)/beta < epsilon) {
				ed50 <- (5-alpha1)/beta1	# ED50
				se <- sqrt(1/ww+(ed50-xbar)^2/ssnwx)/abs(beta1)
				g <- (qnorm(0.975)/beta1)^2/ssnwx
				cl <- ed50+g*(ed50-xbar)/(1-g)
				pm = qnorm(0.975)/beta1/(1-g)*sqrt((1-g)/ww+(ed50-xbar)^2/ssnwx)
				return(list(alpha=alpha1, beta=beta1, ED50=ed50, SE=se, cll =cl-pm, clu=cl+pm))
			}
			alpha <- alpha1
			beta <- beta1
		}
		stop("収束しませんでした")
	}

	# 適合度の検定
	GoodnessOfFitness <- function(x, n, r, alpha, beta)
	{
		p9 <- pnorm(alpha+beta*x-5)
		p0 <- r/n
		cllf <- sum(r*log(p9)+(n-r)*log(1-p9))		# あてはめ後の対数尤度
		q <- 1-p9
		chi <- sum(n*(p0-p9)^2/(q-q^2))			# カイ二乗値
		df <- length(x)-2				# 自由度
		p <- pchisq(chi, df, lower.tail=FALSE)		# P 値
		list(tbl=cbind(x, n, r, n*p9, p0, p9), cllf=cllf, chi=chi, df=df, p=p)
	}

	# 書式付きプリント関数
	printf <- function(fmt, ...) cat(sprintf(fmt, ...))

# ed50 関数本体開始

	transform <- match.arg(transform)			# 引数の補完
	dose <- "x"						# 用量のラベル
	inv <- function(x) x					# 無変換のときの逆関数（無駄だけど）
	if (transform != "none") {
		stopifnot(x > 0)
		if (transform == "ln") {
			x <- log(x)				# 自然対数変換
			dose <- "ln(x)"				# 用量のラベル
			inv <- function(x) exp(x)		# 自然対数の逆関数
		}
		else if (transform == "log") {
			x <- log10(x)				# 常用対数変換
			dose <- "log10(x)"			# 用量のラベル
			inv <- function(x) 10^x			# 常用対数の逆関数
		}
	}

	result <- estimation1(x, n, r)				# 第一近似解
	ollf <- result$ollf					# 観察値の対数尤度
	alpha <- result$alpha					# 切片
	beta <- result$beta					# 傾き

	result <- estimation2(x, n, r, alpha, beta,		# 第二近似解
				epsilon=1e-13)
	alpha <- result$alpha					# 切片
	beta <- result$beta					# 傾き
	printf("P_hat = %g + %g * %s\n", alpha, beta, dose)
	printf("ED50 = %g\n", inv(result$ED50))
	printf("95%% 信頼区間 = [ %g, %g ]\n\n", inv(result$cll), inv(result$clu))

	result <- GoodnessOfFitness(x, n, r, alpha, beta)
	printf("%7s %7s %7s %10s %10s %10s\n", "用量", "n", "r", "e", "r/n", "e/n")
	for (i in 1:nrow(result$tbl)) {
		temp <- result$tbl[i,]
		printf("%7.5g %7.0f %7.0f %10.5g %10.5f %10.5f\n", temp[1], temp[2], temp[3], temp[4], temp[5], temp[6])
	}
	printf("\n")
	printf("カイ二乗値 = %g, 自由度 = %d, P 値 = %.3f\n\n", result$chi, as.integer(result$df), result$p)
	printf("　　観察値の対数尤度 = %g\n", ollf)
	printf("あてはめ後の対数尤度 = %g\n", result$cllf)
}
# 対照群と実験群のケース数，平均値，標準偏差から effect size を求め，統合した effect size を求める
effect.size <- function( Ne, Me, SDe,					# 実験群のケース数，平均値，標準偏差
			 Nc, Mc, SDc)					# 対照群のケース数，平均値，標準偏差
{
	N <- Ne+Nc							# 合計例数
	SD <- sqrt(((Ne-1)*SDe^2+(Nc-1)*SDc^2)/(N-2))			# プールした標準偏差
	g <- (Me-Mc)/SD							# Glass の effect size g
	d <- (1-3/(4*N-9))*g						# g の不偏推定値
	SSest <- N/(Ne*Nc)+d^2/(2*N)					# d の分散
	DF <- data.frame(Ne, Me, SDe, Nc, Mc, SDc, SD, g, d, SSest)	# 結果をデータフレームとしてまとめる
	LST <- list(DF=DF,
		    "mean(g)"=mean(g),					# 統合した effect size（単純平均）
		    "mean(d)"=mean(d),					# 同上（単純平均）
		    "mean(d+)"=sum(d/SSest)/sum(1/SSest))		# 同上（サンプルサイズで重み付け平均）
	class(LST) <- c("effect.size", "list")
	return(LST)
}
# print メソッド（結果のデータフレームを LaTeX 形式で出力し，統合結果も出力する）
print.effect.size <- function(ans)
{
	print.latex(ans[[1]], ctable=FALSE, format="s i i i i i i f6 f6 f6 f6")
	cat(sprintf("g = %.6f\nd = %.6f\nd+ = %.6f\n", ans[[2]], ans[[3]], ans[[4]]))
}
# 標本相関係数の同等性の検定
eq.cor <- function(	n,				# 標本サイズのベクトル
			r)				# 標本相関係数のベクトル
{
	stopifnot(	n > 3,				# 各群の標本サイズは 3 以上でなければならない
			length(n) == length(r),		# n と r の要素数は同じでなければならない
			floor(n) == n,			# 標本サイズは整数値でなくてはならない
			abs(r) <= 1)			# 標本相関係数は -1 ～ 1 の範囲でなければならない
	method <- "標本相関係数の同等性の検定"
	data.name <- paste(deparse(substitute(n)), "and", deparse(substitute(r)))
	k <- length(n)					# 標本の個数
	v <- n-3
	z <- atanh(r)
	sv <- sum(v)
	svz <- sum(v*z)
	chi <- sum(v*z*z)-svz^2/sv			# カイ二乗分布に従う検定統計量
	df <- k-1
	p <- pchisq(chi, df, lower.tail=FALSE)		# P 値	
	result <- list(statistic=c("chi sq."=chi), parameter=c(df=df), p.value=p, method=method, data.name=data.name)
	if (p > 0.05) {					# 標本相関係数が同等と見なせる場合は，点推定値を計算する
		est <- tanh(svz/sv)
		result <- c(result, list(estimate=c("Estimated rho"=est)))
	}
	return(structure(result, class="htest"))
}
# 分散・共分散行列の同等性の検定
eq.cov <- function(	x,			# 第 1 群のデータ行列
			y)			# 第 2 群のデータ行列
{
	x <- subset(x, complete.cases(x))	# 欠損値を持つケースを除く
	y <- subset(y, complete.cases(y))	# 欠損値を持つケースを除く
	p <- ncol(x)				# 変数の個数
	s1 <- var(x)				# 第 1 群の分散・共分散行列
	s2 <- var(y)				# 第 2 群の分散・共分散行列
	n1 <- nrow(x)-1				# 第 1 群の分散・共分散の自由度
	n2 <- nrow(y)-1				# 第 2 群の分散・共分散の自由度
	sa <- (n1*s1+n2*s2)/(n1+n2)		# プールした分散・共分散行列
	chi <- (1-(1/n1+1/n2-1/(n1+n2))*(2*p^2+3*p-1)/(6*p+6))*
		((n1+n2)*log(det(sa))-n1*log(det(s1))-n2*log(det(s2)))
	df <- p*(p+1)/2				# 自由度
	P <- pchisq(chi, df, lower.tail=FALSE)	# P 値
	return(c("Statistics"=chi, "d.f."=df, "P value"=P))
}
euclid <- function(m, n)		# 2^52 未満の 2 つの整数
{
	limit <- 2^52
	stopifnot(is.numeric(m) && is.numeric(n) && m == floor(m) && n == floor(n) && m < limit && n < limit)
	m0 <- m <- abs(m)
	n0 <- n <- abs(n)
	while ((temp <- n %% m) != 0) {
		n <- m
		m <- temp
	}
	lcm <- (m0/m)*n0
	return(list(GCM=m, quotient=c(m0/m, n0/m), LCM=ifelse(lcm > limit, NA, lcm)))
}
eucrid <- euclid
# クラスカル・ウォリス検定（exact test）
exact.kw <- function(	x,							# 分割表（合計を含まない） もしくはデータベクトル
			y=NULL,							# x がデータベクトルのときは，factor ベクトル
			exact=TRUE,						# 正確検定を行うかどうか
			hybrid=FALSE,						# TRUE にすれば，シミュレーションによる
			loop=10000)						# シミュレーションの回数
{

	found <- function()							# 周辺度数が同じ分割表が見つかった
	{
		hh <- sum((um%*%q)^2/rt)					# kw_test(um)
		if (hh >= stat_val || abs(hh-stat_val) <= EPSILON) {
			nprod <- sum(perm_table[rt+1])-sum(perm_table[um+1])
			nntrue <<- nntrue+exp(nprod-nntrue2*log_expmax)
			while (nntrue >= EXPMAX) {
				nntrue <<- nntrue/EXPMAX
				nntrue2 <<- nntrue2+1
			}
		}
		ntables <<- ntables+1
	}

	search <- function(x, y)						# 分割表の探索
	{
		if (y == 1) {							# 見つかった
			found()
		}
		else if (x == 1) {
			if ((t <- um[1, 1]-um[y, 1]) >= 0) {
				um[1, 1] <<- t
				search(nc, y-1)
				um[1, 1] <<- um[1, 1]+um[y, 1]
			}
		}
		else {
			search(x-1, y)
			while (um[y, 1] && um[1, x]) {
				um[y, x] <<- um[y, x]+1
				um[y, 1] <<- um[y, 1]-1
				um[1, x] <<- um[1, x]-1
				search(x-1, y)
			}
			um[y, 1] <<- um[y, 1]+um[y, x]
			um[1, x] <<- um[1, x]+um[y, x]
			um[y, x] <<- 0
		}
	}

	exact.test <- function()						# 正確検定
	{
		denom2 <- 0
		denom <- perm_table[n+1]-sum(perm_table[ct+1])
		while (denom > log_expmax) {
			denom <- denom-log_expmax
			denom2 <- denom2+1
		}
		denom <- exp(denom)
		um[,1] <<- rt
		um[1,] <<- ct
		search(nc, nr)
		printf("正確な P 値 = %.10g\n", nntrue/denom*EXPMAX^(nntrue2-denom2))
		printf("査察した分割表の個数は %s 個\n", ntables)
	}

	kw.test <- function(u)							# クラスカル・ウォリス検定
	{
		return(sum((u%*%q)^2/rt))
	}
	
	monte.carlo <- function()						# モンテカルロ検定
	{
		printf("%i 回のシミュレーションによる P 値 = %g\n", loop, sum(sapply(r2dtable(loop, rt, ct), kw.test) >= stat_val)/loop)
	}

	asymptotic <- function()
	{
		chisq <- (stat_val*12/(n*(n+1))-3*(n+1))/(1-sum(ct^3-ct)/(n^3-n))
		printf("カイ二乗値 = %g, 自由度 = %i, P 値 = %g\n", chisq, nr-1, pchisq(chisq, nr-1, lower.tail=FALSE))
	}

	if (is.list(x)) {
		y <- factor(rep(1:length(x), sapply(x, length)))
		t <- table(y, unlist(x))
	}
	else if (is.matrix(x)) {
		t <- x
	}
	else {
		t <- table(y, x)
	}

	EPSILON <- 1e-10
	EXPMAX <- 1e100	
	log_expmax <- log(EXPMAX)
	nr <- nrow(t)								# 分割表の行数
	nc <- ncol(t)								# 分割表の列数
	rt <- rowSums(t)							# 分割表の行和
	ct <- colSums(t)							# 分割表の列和
	n <- sum(t)								# 総和
	q <- cumsum(c(0,ct[-nc]))+(ct+1)*0.5
	half <- (n+1)*0.5
	stat_val <- kw.test(t)							# クラスカル・ウォリス検定統計量
	asymptotic()								# 検定結果を出力
	if (exact) {
		if (hybrid) {							# モンテカルロ法による検定
			monte.carlo()
		}
		else {								# 正確な検定
			perm_table <- cumsum(c(0, log(1:(n+1))))
			ntables <- nntrue <- nntrue2 <- 0
			um <- matrix(0, nr, nc)
			exact.test()
		}
	}
}
# マン・ホイットニー検定（exact test）
exact.mw <- function(	x,							# 分割表（合計を含まない） もしくはデータベクトル
			y=NULL,							# x がデータベクトルのときは，データベクトル
			exact=TRUE,						# 正確検定を行うかどうか
			hybrid=FALSE,						# TRUE にすれば，シミュレーションによる
			loop=10000)						# シミュレーションの回数
{

	found <- function()							# 周辺度数が同じ分割表が見つかった
	{
		hh <-abs(temp2-sum(um[1,]*r)) # u.test(um)
		if (hh >= stat_val || abs(hh-stat_val) <= EPSILON) {
			nprod <- sum(perm_table[rt+1])-sum(perm_table[um+1])
			nntrue <<- nntrue+exp(nprod-nntrue2*log_expmax)
			while (nntrue >= EXPMAX) {
				nntrue <<- nntrue/EXPMAX
				nntrue2 <<- nntrue2+1
			}
		}
		ntables <<- ntables+1
	}

	search <- function(x, y)						# 分割表の探索
	{
		if (y == 1) {							# 見つかった
			found()
		}
		else if (x == 1) {
			if ((t <- um[1, 1]-um[y, 1]) >= 0) {
				um[1, 1] <<- t
				search(nc, y-1)
				um[1, 1] <<- um[1, 1]+um[y, 1]
			}
		}
		else {
			search(x-1, y)
			while (um[y, 1] && um[1, x]) {
				um[y, x] <<- um[y, x]+1
				um[y, 1] <<- um[y, 1]-1
				um[1, x] <<- um[1, x]-1
				search(x-1, y)
			}
			um[y, 1] <<- um[y, 1]+um[y, x]
			um[1, x] <<- um[1, x]+um[y, x]
			um[y, x] <<- 0
		}
	}

	exact.test <- function()						# 正確検定
	{
		denom2 <- 0
		denom <- perm_table[n+1]-sum(perm_table[ct+1])
		while (denom > log_expmax) {
			denom <- denom-log_expmax
			denom2 <- denom2+1
		}
		denom <- exp(denom)
		um[,1] <<- rt
		um[1,] <<- ct
		search(nc, nr)
		printf("正確な P 値 = %.10g\n", nntrue/denom*EXPMAX^(nntrue2-denom2))
		printf("査察した分割表の個数は %s 個\n", ntables)
	}

	u.test <- function(t)							# マン・ホイットニーの U 検定
	{
		return(abs(temp2-sum(t[1,]*r)))
	}

	monte.carlo <- function()						# モンテカルロ検定
	{
		printf("%i 回のシミュレーションによる P 値 = %g\n", loop, sum(sapply(r2dtable(loop, rt, ct), u.test) >= stat_val)/loop)
	}

	asymptotic <- function()
	{
		z <- stat_val/sqrt(n1n2/(n*(n-1))*(n^3-n-sum(ct^3-ct))/12)
		printf("U = %g, Z = %g, P 値 = %g\n", n1n2/2-stat_val, z, pnorm(z, lower.tail=FALSE)*2)
	}

	if (is.matrix(x)) {							# 分割表が与えられたとき
		t <- x
	}
	else {									# 2 変数が与えられたとき
		nx <- length(x)
		ny <- length(y)
		t <- table(rep(1:2, c(nx, ny)), c(x, y))
	}
	EPSILON <- 1e-10
	EXPMAX <- 1e100
	log_expmax <- log(EXPMAX)
	nr <- nrow(t)								# 分割表の行数
	stopifnot(nr==2)							# 2×k 分割表でないといけない
	nc <- ncol(t)								# 分割表の列数
	rt <- rowSums(t)							# 行和
	ct <- colSums(t)							# 列和
	n1 <- rt[1]								# 第 1 群のケース数
	n2 <- rt[2]								# 第 2 群のケース数
	n1n2 <- n1*n2
	n <- n1+n2								# 全ケース数
	r <- cumsum(c(0,ct[-nc]))+(ct+1)/2
	temp <- n1n2+n1*(n1+1)/2
	temp2 <- temp-n1n2/2
	stat_val <- abs(temp2-sum(t[1,]*r))					# 観察された分割表に対する統計量
	asymptotic()								# 検定結果を出力
	if (exact) {
		if (hybrid) {							# モンテカルロ法による検定
			monte.carlo()
		}
		else {								# 正確な検定
			perm_table <- cumsum(c(0, log(1:(n+1))))
			ntables <- nntrue <- nntrue2 <- 0
			um <- matrix(0, nr, nc)
			exact.test()
		}
	}
}
# 一元配置分散分析（並べ替え検定）
exact.oneway.test <- function(	x,						# リスト
				permutation=TRUE,				# 並べ替え検定を行うかどうか
				hybrid=FALSE,					# TRUE にすれば，シミュレーションによる
				loop=10000)					# シミュレーションの回数
{

	printf <- function(fmt, ...)						# 書式付き出力
	{
		cat(sprintf(fmt, ...))
	}

	found <- function()							# 周辺度数が同じ分割表が見つかった
	{
		hh <- perform.test(um)
		if (hh <= p.value+EPSILON) {
			nprod <- sum(perm_table[rt+1])-sum(perm_table[um+1])
			nntrue <<- nntrue+exp(nprod-nntrue2*log_expmax)
			while (nntrue >= EXPMAX) {
				nntrue <<- nntrue/EXPMAX
				nntrue2 <<- nntrue2+1
			}
		}
		ntables <<- ntables+1
	}

	search <- function(x, y)						# 分割表の探索
	{
		if (y == 1) {							# 見つかった
			found()
		}
		else if (x == 1) {
			if ((t <- um[1, 1]-um[y, 1]) >= 0) {
				um[1, 1] <<- t
				search(nc, y-1)
				um[1, 1] <<- um[1, 1]+um[y, 1]
			}
		}
		else {
			search(x-1, y)
			while (um[y, 1] && um[1, x]) {
				um[y, x] <<- um[y, x]+1
				um[y, 1] <<- um[y, 1]-1
				um[1, x] <<- um[1, x]-1
				search(x-1, y)
			}
			um[y, 1] <<- um[y, 1]+um[y, x]
			um[1, x] <<- um[1, x]+um[y, x]
			um[y, x] <<- 0
		}
	}

	permutation.test <- function()						# 並べ替え検定
	{
		denom2 <- 0
		denom <- perm_table[n+1]-sum(perm_table[ct+1])
		while (denom > log_expmax) {
			denom <- denom-log_expmax
			denom2 <- denom2+1
		}
		denom <- exp(denom)
		um[,1] <<- rt
		um[1,] <<- ct
		search(nc, nr)
		p.value <- nntrue/denom*EXPMAX^(nntrue2-denom2)
		printf("並べ替え検定による P 値 = %.10g\n", p.value)
		printf("査察した分割表の個数は %s 個\n", ntables)
		return(p.value)
	}

	perform.test <- function(u)						# 並べ替え検定
	{
		x <- NULL
		for (i in 1:nr) {
			x <- c(x, rep(score, u[i,]))
		}
		return(oneway.test(x~rep(1:nr, rt))$p.value)
	}
	
	monte.carlo <- function()						# モンテカルロ検定
	{
		p.value <- mean(sapply(r2dtable(loop, rt, ct), perform.test) <= p.value+EPSILON)
		printf("%i 回のシミュレーションによる P 値 = %g\n", loop, p.value)
		return(p.value)
	}

	simple.oneway.test <- function()
	{
		printf("観察値による一元配置分散分析の P 値 = %g\n", p.value)
	}

###### 関数本体

	if (is.list(x)) {
		y <- factor(rep(1:length(x), sapply(x, length)))
		t <- table(y, unlist(x))
	}
	else {
		stop("群ごとのデータはリストで与えること")
	}

	EPSILON <- 1e-10
	EXPMAX <- 1e100	
	log_expmax <- log(EXPMAX)
	nr <- nrow(t)								# 分割表の行数
	nc <- ncol(t)								# 分割表の列数
	rt <- rowSums(t)							# 分割表の行和
	ct <- colSums(t)							# 分割表の列和
	n <- sum(t)								# 総和
	q <- cumsum(c(0,ct[-nc]))+(ct+1)*0.5
	half <- (n+1)*0.5
	score <- as.numeric(colnames(t))
	p.value <- perform.test(t)						# 観察値による検定
	simple.oneway.test()								# 検定結果を出力
	if (permutation) {
		if (hybrid) {							# モンテカルロ法による検定
			p.value <- monte.carlo()
		}
		else {								# 並べ替え検定
			perm_table <- cumsum(c(0, log(1:(n+1))))
			ntables <- nntrue <- nntrue2 <- 0
			um <- matrix(0, nr, nc)
			p.value <- permutation.test()
		}
	}
	invisible(p.value)
}
excel.w <- function(nc)
{
	matrix(scan("clipboard", quiet=TRUE), byrow=TRUE, ncol=nc)
}
excel <- function(nc)
{
	matrix(scan("", quiet=TRUE), byrow=TRUE, ncol=nc)
}
extended.chisq.test <- function (tbl) 
{
	method <- "二次元以上の集計表での独立性の検定"
	data.name <- deparse(substitute(tbl))
    n <- sum(tbl)
    dm <- dim(tbl)
    variables <- length(dm)
    if (variables > 1) {
        m <- vector("list", length=variables)
        for (k in 1:variables) m[[k]] <- apply(tbl, k, sum)/n
        expected <- apply(expand.grid(m), 1, prod)*n
        chi.sq <- sum((c(tbl)-expected)^2/expected)
        df <- prod(dm)-1-sum(dm-1)
        p.value <- pchisq(chi.sq, df, lower.tail=FALSE)
        if (any(expected < 5)) {
        	warning("カイ自乗近似は不正確かもしれません")
        }
        return(structure(list(statistic=c("X-squared"=chi.sq), 
        	parameter=c(df=df), p.value=p.value, method=method,
        	data.name=data.name, tbl=tbl), class="htest"))
    }
    else {
	    cat("二次元以上の集計表が対象です")
	 }
}
# 因子分析の適合度検
fa.fit.test <- function(data,					# データ行列（データフレーム）
			factors,				# 抽出した因子数
			uniquenesses)				# 独自性ベクトル
{
	p <- ncol(data)						# 変数の個数
	n <- nrow(data)						# ケース数
	r <- cor(data)						# 相関係数行列
	sc <- diag(1/sqrt(uniquenesses))			# 
	r <- sc%*%r%*%sc					# 
	e <- eigen(r, symmetric=TRUE, only.values=TRUE)		# 固有値だけ求める
	e <- e$values[-(1:factors)]				# 抽出された因子で説明されない部分
	s <- -sum(log(e)-e)+factors-p				# 統計量
	chisq <- (n-1-(2*p+5)/6-(2*factors)/3)*s		# カイ二乗統計量に変換
	df <- ((p-factors)^2-p-factors)/2			# 自由度
	P <- pchisq(chisq, df, lower.tail=FALSE)		# P 値
	return(c(Chi.sq.=chisq, df=df, P=P))			# 名前付ベクトルで返す
}
# チャーノフの顔グラフで使用するデータを調整する
face.data <- function(	d,								# データ行列
			pos=1:18)							# 18 個のパーツに対応させるデータ行列の列番号（変数）
{
	lo <- c(rep(0.2, 5), 0.1, 0.2, -5, 0.2, 0.1, 0.1, 0.3, 0.1, 0.3, rep(0.1, 4))
	hi <- c(0.8, 0.8, 1, 0.8, 0.8, 0.4, 0.8, 5, 0.8, 0.7, 0.9, 0.7, rep(0.9, 4), 1, 0.9)
	fx <- c(0.5, 0.5, 1, 0.5, 0.5, 0.2, 0.5, 0, 0.5, 0.4, rep(0.5, 3), 0.6, 0.5, 0.5, 1, 0.5)
	n <- nrow(d)									# ケース数
	min.d <- apply(d, 2, min)
	max.d <- apply(d, 2, max)
	x <- matrix(0, n, 18)
	for (i in 1:18) {
		k <- pos[i]
		if (k > 0) {
			min.dk <- min.d[k]
			max.dk <- max.d[k]
			if (min.dk == max.dk) {						# 全部が同じ値を持つようなとき，
				x[,i]  <- fx[i]						# 平均的なパラメータ
			}
			else {
				x[,i] <- (d[,k]-min.dk)/(max.dk-min.dk)*(hi[i]-lo[i])+lo[i]
			}
		}
	}
	return(x)
}
# チャーノフの顔グラフ
face.plot <- function(x, size=480)
{
	# 円弧を描く
	arc1 <- function(x1, y1, r, l)
	{
		sign <- ifelse(l > 0, -1, 1)
		theta <- sign*acos(x1/r)
		y1 <- y1-sign*sqrt(r^2-x1^2)
		if (l <= 0) {
			arc(0, y1, r, theta, pi-theta)
		}
		else {
			arc(0, y1, r, pi-theta, pi*2+theta)
		}
	}
	# 円（の一部）を描く
	arc <- function(ox, oy, r, theta.start, theta.end)
	{
		step <- min(0.1, (theta.end-theta.start)*0.1)
		interval <- c(seq(theta.start, theta.end, step), theta.end)
		lines(r*cos(interval)+ox, r*sin(interval)+oy)
	}
	# 楕円（の一部）を描く
	ellipse <- function(ox, oy, r.a, r.b, theta.axis, theta.start, theta.end)
	{
		theta.end <- theta.end+(theta.end <= theta.start)*pi*2
		temp1 <- r.a*r.b
		temp2 <- 30/(r.a+r.b)
		k <- (theta.end-theta.start)/temp2+2
		x <- y <- numeric(k)
		for (i in 1:(k-1)) {
			factor <- temp1/sqrt((r.a*sin(theta.start))^2+(r.b*cos(theta.start))^2)
			x[i] <- factor*cos(theta.axis+theta.start)
			y[i] <- factor*sin(theta.axis+theta.start)
			theta.start <- theta.start+temp2
		}
		factor <- temp1/sqrt((r.a*sin(theta.end))^2+(r.b*cos(theta.end))^2)
		x[k] <- factor*cos(theta.axis+theta.end)
		y[k] <- factor*sin(theta.axis+theta.end)
		lines(ox+x, oy+y)
	}

	pi2 <- 2*pi
	plot(c(-500, 500), c(-500, 500), type="n", xlab="", xaxt="n", ylab="", yaxt="n", bty="n")
	size2 <- size*(1+x[1])/2
	theta <- (pi/4)*(2*x[2]-1)
	h <- size*(1+x[3])/2
	x1 <- size2*cos(theta)
	y1 <- size2*sin(theta)
# 顔の上半分 
	ak <- 1-x[4]^2
	oy1 <- (ak*x1^2+y1^2-h^2)/(2*(y1-h))
	r.b1 <- h-oy1
	r.a1 <- r.b1/sqrt(ak)
	theta.start <- atan((y1-oy1)/x1)
	theta.end <- pi-theta.start
	ellipse(0, oy1, r.a1, r.b1, 0, theta.start, theta.end)
# 顔の下半分 
	ak <- 1-x[5]^2
	oy2 <- (ak*x1^2+y1^2-h^2)/(2*(y1+h))
	r.b2 <- h+oy2
	r.a2 <- r.b2/sqrt(ak)
	theta.end <- atan((y1-oy2)/x1)
	theta.start <- pi-theta.end
	ellipse(0, oy2, r.a2, r.b2, 0, theta.start, theta.end)
# 鼻 
	y <- h*x[6]
	lines(c(0, 0), c(y, -y))
# 口 
	pm <- -h*(x[7]+(1-x[7])*x[6])
	wm <- sqrt(r.a2^2*(1-(pm-oy2)^2/r.b2^2))
	if (x[8] == 0) {
		lines(c(-wm/2, wm/2), c(pm, pm))
	}
	else {
		r <- h/abs(x[8])
		am <- x[9]*r
		x1 <- ifelse(am > wm, x[9]*wm, am)
		l <- ifelse(x[8] <= 0, -1, 1)
		y1 <- pm-l*(r-sqrt(r^2-x1^2))
		arc1(x1, y1, r, l)
	}
# 目 
	ye <- h*(x[10]+(1-x[10])*x[6])
	we <- sqrt(r.a1^2*(1-(ye-oy1)^2/r.b1^2))
	xe <- we*(1+2*x[11])/4
	theta <- (2*x[12]-1)*pi/5
	r.a3 <- x[14]*min(xe, we-xe)
	r.b3 <- sqrt(r.a3^2*(1-x[13]^2))
	ellipse(xe, ye, r.a3, r.b3, theta, 0, pi2)
	ellipse(-xe, ye, r.a3, r.b3, pi-theta, 0, pi2)
# 瞳 
	re <- r.a3/sqrt(cos(theta)^2+sin(theta)^2/x[13]^2)
	shift <- re*(2*x[15]-1)
	sapply(c(xe, -xe)-shift, function(arg) arc(arg, ye, 3, 0, pi2))
# 眉 
	theta2 <- 2*(1-x[17])*(pi/5)
	theta3 <- ifelse(theta >= 0, theta+theta2, theta-theta2)
	len <- re*(2*x[18]+1)/2
	x0 <- len*cos(theta3)
	x1 <- xe-c(x0, -x0)
	y0 <- len*sin(theta3)
	y1 <- ye+2*(x[16]+0.3)*r.a3*x[13]-c(y0, -y0)
	lines(x1-shift, y1)
	lines(-x1-shift, y1)
}
# factanal 関数の結果を整形して表示する
factanal2 <- function(	dat,							# データ行列
			factors=0,						# 抽出する因子数
			rotation=c("promax", "varimax", "none"),		# 因子軸の回転法
			scores=c("none", "regression", "Bartlett"),		# 因子得点の算出法
			verbose=TRUE)						# 結果の画面表示をするかどうか
{
	p <- ncol(dat)								# 変数の個数
	n <- nrow(dat)								# 行数（欠損値を含むケースも含む）
	if (is.null(colnames(dat))) colnames(dat) <- paste("Var", 1:p, sep=".")	# 変数名がないときにデフォルト名をつける
	if (is.null(rownames(dat))) rownames(dat) <- paste("Case", 1:n, sep="-")# 行名がないときにデフォルト名をつける
	dat <- subset(dat, complete.cases(dat))					# 欠損値を持つケースを除く
	rotation <- match.arg(rotation)						# 引数の補完
	scores <- match.arg(scores)						# 引数の補完
	
	if (factors == 0) {							# 抽出因子数が指定されないときは，
		factors <- max(1, floor((2*p+1-sqrt((2*p+1)^2-4*(p^2-p)))/2))	# デフォルトの因子数
	}
	txt <- sprintf('factanal(dat, factors=factors, rotation="%s", scores=scores)', rotation) # rotation に実際の文字列を渡すためにこのようにする
	result <- eval(parse(text=txt))						# 関数呼び出し
	Communality <- 1-result$uniquenesses					# 共通性は，斜交回転のときには因子負荷量の二乗和ではない
	result$cosmetic <- cbind(result$loadings, Communality)			# 共通性を付加
	if (rotation!="promax") {						# 斜交回転でない場合には，
		SS.loadings <- colSums(result$loadings^2)			# 因子負荷量の二乗和
		SS.loadings <- c(SS.loadings, sum(SS.loadings))			# 総和を加える
		Proportion <- SS.loadings/p*100					# 寄与率
		Cum.Prop. <- cumsum(Proportion)					# 累積寄与率
		Cum.Prop.[factors+1] <- NA
		result$cosmetic <- rbind(result$cosmetic, SS.loadings, Proportion, Cum.Prop.)
	}
	if (verbose == TRUE) {							# 画面表示をするとき
		if (result$dof) {						# モデル適合度の自由度が 0 でないとき
			test <- data.frame(result$STATISTIC, result$dof, result$PVAL)
			colnames(test) <- c("Chi sq.", "d.f.", "P value")
			rownames(test) <- ""
			cat(sprintf("H0: %i factors are sufficient.\n", factors))
			print(test)
		}
		else {								# 自由度が 0 になるとき
			cat(sprintf("The degrees of freedom for the model is 0 and the fit was %g\n", result$criteria[1]))
		}
		cat(sprintf("\nFactor loadings(rotation:%s)\n", rotation))	# 因子負荷量
		print(result$cosmetic)
		if (scores != "none") {
			cat(sprintf("\nFactor scores(%s)\n", scores))		# 因子得点
			print(result$scores)
		}
	}
	invisible(result)							# 明示的に print しないと，何も出力しない
}
# プロマックス解の因子間相関係数行列
factor.correlation <- function(x, factors, ...)
{
    ans <- factanal(x, factors, rotation="none", ...)		# 回転を行わない結果を求める
    ans2 <- promax(ans$loadings)				# プロマックス回転による結果を求める
    name <- colnames(ans2$loadings)				# 名前の保存
    o <- order(colSums(ans2$loadings^2), decreasing=TRUE)	# SS loadings の大きい順
    ans2$loadings <- ans2$loadings[, o]				# loadings の並べ替え（行）
    colnames(ans2$loadings) <- name				# 名前の付け替え
    class(ans2$loadings) <- "loadings"				# class がなくなるので再設定
    ans2$rotmat <- ans2$rotmat[o, o]				# rotmat の並べ替え（行・列）
    ans3 <- ans2$rotmat						# 回転行列を取り出す
    r <- solve(t(ans3) %*% ans3)				# 因子間相関係数行列を計算する
    colnames(r) <- rownames(r) <- name				# 名前を付ける（必須ではない）
    return(list(loadings=ans2$loadings, r=r))			# プロマックス解と因子間相関係数行列
}
# 多重共線性のある変数を検出する
# http://www.ec.kagawa-u.ac.jp/~hori/spss/eig0vec.txt
find.multico <- function(	data,				# データ行列
				epsilon=1e-10)			# 限界値
{
	result <- eigen(cor(data[complete.cases(data),]))	# 固有値・固有ベクトルを求める
	n <- length(result$values)
	result <- matrix(c("", "***")[(abs(result$vectors[,result$values < epsilon]) > epsilon)+1], nr=n)
	if (ncol(result)) {
		rownames(result) <- colnames(data)
		result
	}
	else {
		"no multico"
	}
}
# 有限母集団の割合の推定に必要な標本サイズを決める
finite <- function(	n=NULL,			# 標本サイズ
			N=NULL,			# 母集団サイズ
			p=NULL,			# 母比率
			epsilon=NULL,		# 精度
			conf.level=NULL)	# 信頼度（信頼係数）
{
	if (sum(sapply(list(n, N, p, epsilon, conf.level), is.null)) != 1) {
		stop("n, N, p, epsilon, conf.level のうちの，どれか一つだけが NULLでなければならない")
	}
	n.function <- quote(N/((epsilon/qnorm(0.5-conf.level/2, lower.tail=FALSE))^2*((N-1)/(p*(1-p)))+1))
	if (is.null(n)) {
		n <- eval(n.function)
	}
	else if (is.null(epsilon)) {
		epsilon <- uniroot(function(epsilon) eval(n.function)-n, c(1e-7, 0.9999999))$root
	}
	else if (is.null(N)) {
		N <- uniroot(function(N) eval(n.function)-n, c(1, 1e7))$root
	}
	else if (is.null(p)) {
		p <- uniroot(function(p) eval(n.function)-n, c(1e-7, 0.9999999))$root
	}
	else if (is.null(conf.level)) {
		conf.level <- uniroot(function(conf.level) eval(n.function)-n, c(1e-7, 0.9999999))$root
	}
	METHOD <- "有限母集団の割合の推定に必要な標本サイズ"
	structure(list(n=n, N=N, p=p, epsilon=epsilon, conf.level=conf.level, method=METHOD), class="power.htest")
}
# Fisher の exact test
fisher <- function(	x,							# 分割表（合計を含まない） もしくは factor ベクトル
			y=NULL,							# x が factor ベクトルのときは，factor ベクトル
			exact=TRUE,						# 正確検定を行うかどうか
			method=c("Fisher", "Pearson"),				# フィッシャーによるか，ピアソンによるか
			hybrid=FALSE,						# TRUE にすれば，シミュレーションによる
			loop=10000)						# シミュレーションの回数
{

	found <- function()							# 周辺度数が同じ分割表が見つかった
	{
		if (method == "Fisher") {					# フィッシャーの基準による
			nprod <- temp-sum(perm_table[um+1])
			if (nprod <= criterion+EPSILON) {
				nntrue <<- nntrue+exp(nprod-nntrue2*log_expmax)
				while (nntrue >= EXPMAX) {
					nntrue <<- nntrue/EXPMAX
					nntrue2 <<- nntrue2+1
				}
			}
		}
		else {								# ピアソンの基準による 
			hh <- sum((um-ex)^2/ex) # chisq.test(um)
			if (hh >= stat_val || abs(hh-stat_val) <= EPSILON) {
				nprod <- temp-sum(perm_table[um+1])
				nntrue <<- nntrue+exp(nprod-nntrue2*log_expmax)
				while (nntrue >= EXPMAX) {
					nntrue <<- nntrue/EXPMAX
					nntrue2 <<- nntrue2+1
				}
			}
		}
		ntables <<- ntables+1
	}

	search <- function(x, y)						# 分割表の探索
	{
		if (y == 1) {							# 見つかった
			found()
		}
		else if (x == 1) {
			if ((t <- um[1, 1]-um[y, 1]) >= 0) {
				um[1, 1] <<- t
				search(nc, y-1)
				um[1, 1] <<- um[1, 1]+um[y, 1]
			}
		}
		else {
			search(x-1, y)
			while (um[y, 1] && um[1, x]) {
				um[y, x] <<- um[y, x]+1
				um[y, 1] <<- um[y, 1]-1
				um[1, x] <<- um[1, x]-1
				search(x-1, y)
			}
			um[y, 1] <<- um[y, 1]+um[y, x]
			um[1, x] <<- um[1, x]+um[y, x]
			um[y, x] <<- 0
		}
	}

	exact.test <- function()						# 正確検定
	{
		denom2 <- 0
		denom <- perm_table[n+1]-sum(perm_table[ct+1])
		while (denom > log_expmax) {
			denom <- denom-log_expmax
			denom2 <- denom2+1
		}
		denom <- exp(denom)
		um[,1] <<- rt
		um[1,] <<- ct
		search(nc, nr)
		printf("%s の方法による，正確な P 値 = %.10g\n", method, nntrue/denom*EXPMAX^(nntrue2-denom2))
		printf("査察した分割表の個数は %s 個\n", ntables)
	}

	chisq.test <- function(t)						# カイ二乗近似検定
	{
		return(sum((t-ex)^2/ex))
	}

	prob <- function(t)
	{
		return(temp-sum(perm_table[t+1]))
	}

	monte.carlo <- function()						# モンテカルロ検定
	{
		if (method == "Fisher") {					# フィッシャーの基準による
			count <- sum(sapply(r2dtable(loop, rt, ct), prob) <= criterion+EPSILON)
		}
		else {								# ピアソンの基準による
			count <- sum(sapply(r2dtable(loop, rt, ct), chisq.test) >= stat_val)
		}
		printf("%i 回のシミュレーション（%s の方法）による P 値 = %g\n", loop, method, count/loop)
	}

	asymptotic <- function()
	{
		printf("カイ二乗値 = %g, 自由度 = %i, P 値 = %g\n", stat_val, df, pchisq(stat_val, df, lower.tail=FALSE))
	}

	if (is.matrix(x)) {							# 分割表が与えられたとき
		t <- x
	}
	else {									# 2 変数が与えられたとき
		t <- table(y, x)
	}

	EPSILON <- 1e-10
	EXPMAX <- 1e100
	log_expmax <- log(EXPMAX)
	nr <- nrow(t)								# 分割表の行数
	nc <- ncol(t)								# 分割表の列数
	rt <- rowSums(t)							# 分割表の行和
	ct <- colSums(t)							# 分割表の列和
	n <- sum(t)								# 総和
	ex <- outer(rt, ct)/n							# 期待値
	stat_val <- chisq.test(t)						# 観察された分割表のカイ二乗値
	df <- (nr-1)*(nc-1)							# 自由度
	asymptotic()								# 検定結果を出力
	if (exact) {								# 正確検定を行うなら，
		method <- match.arg(method)					# ピアソンの基準かフィッシャーの基準か
		perm_table <- cumsum(c(0, log(1:(n+1))))			# 対数を取った階乗の表
		temp <- sum(perm_table[rt+1])
		criterion <- temp-sum(perm_table[t+1])
		if (hybrid) {							# モンテカルロ法による検定
			monte.carlo()
		}
		else {								# 正確な検定
			ntables <- nntrue <- nntrue2 <- 0
			um <- matrix(0, nr, nc)
			exact.test()
		}
	}
}
# 度数分布表を与えて，正規分布へのあてはめを行い，指定によっては図を描く
fit.normal <- function(	f,						# 度数を表すベクトル
			l,						# 最小の階級の下限値
			w,						# 階級幅
			accuracy=0,					# 測定精度
			method=c("density", "area"),			# あてはめ方法。確率密度による場合は "density"，面積を計算する場合は "area"
			xlab="x",					# 結果のグラフ表示における横軸の名称
			ylab="f(x)",					# 結果のグラフ表示における縦軸の名称
			plot=TRUE,					# 結果のグラフを描くかどうか
			col="gray",					# ヒストグラムの塗りつぶし色
			col1="blue",					# 理論正規分布曲線の色
			col2="red")					# 期待値を描く○の色
{
	method <- match.arg(method)					# 引数の略語を補完する
	f <- c(0, f, 0)							# 度数ベクトルの前後に階級を一つ付加する
	m <- length(f)							# 階級数
	x <- l-accuracy/2-3*w/2+w*(1:m)					# 級限界
	names(f) <- as.character(x)					# 階級の名称
	n <- sum(f)							# ケース数
	mean <- sum(f*x)/n						# 平均値の推定値
	sd <- sqrt(sum(f*x^2)/n-mean^2)					# 標準偏差の推定値
	if (method == "area") {						# 面積を計算する方法の場合
		z <- (x+w/2-mean)/sd					# 級限界の標準得点
		F <- pnorm(z)						# 下側確率
		F[m] <- 1						# 全部加えたら 1
		p <- F-c(0, F[-m])					# 各階級の確率
	}
	else {								# 確率密度による場合
		z <- (x-mean)/sd					# 標準得点
		p <- dnorm(z)*w/sd					# 確率密度
	}
	exp <- n*p
	if (plot) {							# 結果を図示するときは
		xl <- l+w*-1:(m-1)					# 横軸の値
		plot(xl, c(f, 0)/n, type="n", xlab=xlab, ylab=ylab)	# 図の枠組み
		rect(xl, 0, xl+w, c(f, 0)/n, col=col)			# ヒストグラムの長方形を描く
		points(x, p, col=col2)					# 期待値の点を打つ
		x2 <- seq(x[1], x[m], length=100)			# 理論値
		lines(x2, dnorm(x2, mean, sd)*w, col=col1)		# 理論値の度数多角形を描く
		abline(h=0)						# 横軸を描く
	}
	result <- list(method=method, mean=mean, sd=sd, table=cbind(x, f, z, p, exp))
	class(result) <- "fit.normal"					# 結果にクラス名を与える
	return(result)
}
print.fit.normal <- function(x)						# fit.normal クラスの出力関数
{
	cat("正規分布へのあてはめ (方法 =", x$method, ")\n\n")
	cat(sprintf("推定された　　平均値 = %.7g  標準偏差 = %.7g\n\n", x$mean, x$sd))
	x <- x$table
	cat("   級中心    頻度      z 値       密度      期待値\n")
	for (i in 1:nrow(x)) {
		cat(sprintf("%8.4f %6i %9.3f %9.5f %10.2f\n", x[i,1], as.integer(x[i,2]), x[i,3],  x[i,4], x[i,5]))
	}
	cat(sprintf("     合計 %6i %19.5f %10.2f\n", sum(x[,2]), sum(x[,4]), sum(x[,5])))
}
# コピー，カットによりクリップボードに取り込んだ内容を scan により R に取り込む
fixed <- function()
{
	x <- scan("", sep="\n", what="", quiet=TRUE, strip.white=TRUE)		# 一行ずつ，前後の空白をはぎ取った文字列を読む
	n <- length(x)								# 何行読んだか
	res <- sapply(x, function(x) strsplit(x, "([\t]| +)"))			# 空白類を区切り文字として切り取る
	len <- sapply(res, length)						# 何項目読んだか
	if (len[1] == len[2]-1) {						# 一行目と二行目の項目数が違うとき，
		res[[1]] <- append("", res[[1]])				# 一行目の最初の項目として空を挿入する
		len[1] <- len[1]+1						# 1 項目増える
	}
	max.len <- max(len)							# 最大項目数
	res2 <- matrix("", n, max.len)						# 行列にする
	for (i in 1:n) {
		res2[i,] <- c(unlist(res[i]), rep("", max.len-len[i]))		# 一行の残りの項目は空とする
	}
	class(res2) <- c("fixed", "matrix")					# fixed クラスにする（print.fixed を使うために）
	return(res2)
}
# 度数分布表を作成する
freq <- function(	x,					# データベクトル
			lo,					# 階級値下限
			hi,					# 階級値上限
			w)					# 階級幅
{
	x <- x[!is.na(x)]					# 欠損値を持つケースを除く
	n <- length(x)						# サンプルサイズ
	f <- table(cut(x, br = seq(lo, hi, w), right = FALSE))	# 度数分布を得る
	if (n != sum(f)) {
		stop("lo, hi の指定を変えて，範囲を広げてください。")
	}
	res <- cbind(f, f/n*100, cumsum(f)/n*100)		# 度数，相対度数，累積相対度数
	colnames(res) <- c("度数", "相対度数", "累積相対度数")
	return(res)
}
######
#
# 分析対象変数の度数分布表を作成し，変数が factor のときには棒グラフ，数値変数の場合にはヒストグラムを描く
#
######

frequency <- function(i,							# 分析対象変数のあるデータフレームの列番号または変数名ベクトル
                      df,							# データフレーム
                      latex=TRUE,						# LaTeX 形式で度数分布表を出力する（デフォルトは LaTeX 形式）
                      captions=NULL,						# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
                      labels=NULL,						# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
                      output="",						# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
		      encoding=getOption("encoding"),				# ファイルに出力するときのエンコーディング（デフォルトは OS による）
                      plot="",							# 棒グラフ・ヒストグラムを描き出すファイル名（デフォルトは Quarts デバイスに出力）
                      type=c("pdf", "png", "jpeg", "bmp", "tiff"),		# 画像フォーマット（plot と併せてファイル名の拡張子として使う）
                      width=500,						# 画像の横幅のピクセル数（デフォルトは500ピクセル）
                      height=375,						# 画像の高さのピクセル数（デフォルトは375ピクセル）
                      xlab=NULL,						# 横軸のラベル（デフォルトは対象変数名）。何も描かないときは空文字列を指定する
                      ylab="度数",						# 縦軸のラベル（デフォルトは「度数」）。何も描かないときは空文字列を指定する
                      main="")							# グラフのメインタイトル（デフォルトはメインタイトルを付けない）
{

	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

										# グラフをファイルに出力する場合の処理
	if (plot != "") {							# グラフをファイルに出力するとき plot はファイル名（拡張子を除く）
		type <- match.arg(type)						# 画像ファイルの形式
		if (type == "pdf") {						# pdf なら，一つ一つの画像を別々のファイルに出力するために onefile = FALSE にする
			pdf(sprintf("%s%%03i.pdf", plot), onefile=FALSE, width=width/72, height=height/72)
			   							# pdf は，画像の大きさの指定がインチ単位なので 72dot/inch で換算
		}
		else if (type == "png") {
			png(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else if (type == "bmp") {
			bmp(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else if (type == "jpeg") {
			jpeg(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else { # type == "tiff"
			tiff(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
	}

	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
# それぞれの分析対象変数について
	index <- 0
	for (i1 in i) {								# i には，分析対象とする変数のデータフレーム上での列位置を示す番号がベクトルとして入っている
		x <- subset(df[,i1], !is.na(df[,i1]))				# 欠損値を除いたデータを対象とする
		index <- index+1
		v.name <- colnames(df)[i1]					# 変数の名前を取り出す
		xlab2 = if (is.null(xlab)) v.name else xlab			# x 軸のラベル。デフォルト（NULL）なら集計対象とする変数名。描かないなら空文字列
		
# 変数が factor のとき

		if (is.factor(x)) {						# 集計対象変数が factor なら，度数と相対度数だけを求め，棒グラフ（barplot）を描く
			count <- table(x)					# table で，度数分布を求める
			n <- sum(count)						# NA を除く有効ケース数
			freq <- count/n*100					# 相対度数（%）を求める
			name <- names(freq)					# 各カテゴリーの名前
			ln <- length(freq)					# カテゴリー数
			if (latex) {								# LaTeX 形式で度数分布表を出力する
				cat("\n\\begin{table}[htbp]\n", file=output)			# \begin{table}[htbp]
				if (is.null(captions)) {
					cat(sprintf("\\caption{%sの度数分布}\n", v.name), file=output)	# \caption{○○の度数分布}
				}
				else {
					cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
				}
				if (!is.null(labels)) {
					cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
				}
				cat("\\centering\n", file=output)				# \centering
				cat("\\begin{tabular}{lrr} \\hline\n", file=output)		# \begin{tabular}{lrr} \hline
				cat("項目\t&\t度数\t&\t相対度数\\\\ \\hline\n", file=output)	# 項目 & 度数 & 相対度数 \\ \hline
				for (j in 1:ln) {						# 項目の数だけ繰り返す
					cat(sprintf("%s\t&\t%i\t&\t%.1f\\\\%s\n",		# ○○ & ○○ & ○○.○ \\
					    name[j], count[j], freq[j],				# ○○ & ○○ & ○○.○ \\
					    if (j == ln) "\\hline" else ""), file=output)	# ○○ & ○○ & ○○.○ \\ \hline
				}
				cat(sprintf("合計\t&\t%i\t&\t100.0\\\\ \\hline\n", n),		# 合計 & ○○ & ○○.○ \\ \hline
				    file=output)
				cat("\\end{tabular}\n", file=output)				# \end{tabular}
				cat("\\end{table}\n", file=output)				# \end{table}
			}
			else {									# tab で区切って出力する
				cat("\n表　", v.name, "の度数分布\n\n", sep="", file=output)	# 表　○○の度数分布
				cat("項目", "度数", "相対度数\n", sep="\t", file=output)	# 項目	度数	相対度数
				for (j in 1:ln) {						# 項目の数だけ繰り返す
					cat(name[j], count[j], sprintf("%.1f\n", freq[j]),	# ○○	○○	○○.○
					    sep="\t", file=output)				# ○○	○○	○○.○
				}								# ○○	○○	○○.○
				cat("合計", n, sprintf("%.1f\n", 100), sep="\t", file=output)	# 合計	○○	○○.○
			}
			barplot(count, xlab=xlab2, ylab=ylab)					# 棒グラフ（barplot）を描く
		}

# 数値変数のとき

		else {								# 集計対象変数が factor でない（間隔尺度・比尺度）なら，累積度数も求め，ヒストグラムを描く
			options(warn=-1)					# 次の行のような hist の使い方（plot=FALSE）で，R2.4.0 が過剰な warning を出すので対策を取った（R2.4.1 からは不要）
			ans <- hist(x, right=FALSE, plot=FALSE)			# 階級，度数などの hist オブジェクトを得る。級限界値は，左を含み，右を含まない
			ln <- length(ans$breaks)				# breaks の個数（階級数より 1 だけ少ない）
			if (ans$breaks[ln] == max(x)) {				# もし breaks の一番大きい値が，データ中の最大値なら，その値は右側級限界に含まれてしまっている
				ans <- hist(x, breaks=seq(ans$breaks[1], by=diff(ans$breaks)[1], length=ln+1), right=FALSE, plot=FALSE)
										# 一番大きい階級の上に階級数を一つ増やして度数分布を取り直す
			}
			options(warn=0)						# warning のレベルを元に戻す
			count <- ans$counts					# hist オブジェクトから，度数を取り出す
			ccount <- cumsum(count)					# 累積度数を作る
			n <- sum(count)						# サンプルサイズ
			freq <- count/n*100					# 相対度数（%）
			cfreq <- ccount/n*100					# 累積相対度数（%）
			ln <- length(freq)					# 階級数
			name <- ans$breaks[1:ln]				# 階級の名前は，級限界値とする
			fraction <- 0						# 階級表示に必要な小数点以下の桁数
			for (j in 1:ln) {
				char.vec <- unlist(strsplit(as.character(name[j]), "\\."))
				if (length(char.vec) == 2) {
					fraction <- max(fraction, nchar(char.vec[2]))
				}
			}
			if (fraction == 0) {
				fmt <- "%s～"
			}
			else {
				fmt <- sprintf("%%.%sf～", fraction)
			}
			if (latex) {								# LaTeX 形式で度数分布表を出力する
				cat("\n\\begin{table}[htbp]\n", file=output)			# \begin{table}[htbp]
				if (is.null(captions)) {
					cat(sprintf("\\caption{%sの度数分布}\n", v.name), file=output)	# \caption{○○の度数分布}
				}
				else {
					cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
				}
				if (!is.null(labels)) {
					cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
				}
				cat("\\centering\n", file=output)				# \centering
				cat("\\begin{tabular}{lrrrr} \\hline\n", file=output)		# \begin{tabular}{lrrrr} \hline
				cat("階級\t&\t度数\t&\t相対度数\t&\t累積度数\t&\t累積相対度数\\\\ \\hline\n", file=output)
												# 階級 & 度数 & 相対度数 & 累積度数 & 累積相対度数 \\ \hline
				fmt <- paste(fmt, "\t&\t%i\t&\t%.1f\t&\t%i\t&\t%.1f\\\\%s\n", sep="")
				for (j in 1:ln) {						# 階級の数だけ繰り返す
					cat(sprintf(fmt,
					    name[j], count[j], freq[j], ccount[j], cfreq[j],	# ○○ & ○○ & ○○.○ & ○○ & ○○.○ \\
					    if (j == ln) " \\hline" else ""), file=output)	# ○○ & ○○ & ○○.○ & ○○ & ○○.○ \\ \hline
				}
				cat(sprintf("合計\t&\t%i\t&\t100.0\\\\ \\hline\n", n),		# 合計 & ○○ & ○○.○ & ○○ & ○○.○ \\ \hline
				    file=output)
				cat("\\end{tabular}\n", file=output)				# \end{tabular}
				cat("\\end{table}\n", file=output)				# \end{table}
			}
			else {									# tab で区切って出力する
				cat("\n表　", v.name, "の度数分布\n\n", sep="", file=output)	# 表　○○の度数分布
				cat("階級\t度数\t相対度数\t累積度数\t累積相対度数\n",		# 階級	度数相対度数	累積度数	累積相対度数
				    file=output)
				fmt <- paste(fmt, "\t%i\t%.1f\t%i\t%.1f\n", sep="")
				for (j in 1:ln) {						# 階級の数だけ繰り返す
					cat(sprintf(fmt,					# 階級，度数
					    name[j], count[j], freq[j], ccount[j], cfreq[j]),	# 相対度数，累積度数
					    file=output)					# 累積相対度数
				}
				cat("合計", n, sprintf("%.1f\n", 100), sep="\t", file=output)	# 合計，○○，○○.○，○○，○○.○
			}
			plot(ans, main=main, xlab=xlab2, ylab=ylab)				# ヒストグラムを描く（main, xlab, ylab を指定できる）
		}
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
										# グラフをファイルに出力した場合の後始末
	if (plot != "") {							# ファイルに出力しているなら，
		dev.off()							# デバイスを閉じる
	}
}
# フリードマン検定（plus 多重比較）
friedman <- function(dat)						# データ行列
{
	method <- "フリードマン検定（plus 多重比較）"
	data.name <- deparse(substitute(dat))
	dat <- subset(dat, complete.cases(dat))				# 欠損値を持つケースを除く
	row <- nrow(dat)						# ケース数
	col <- ncol(dat)						# 条件数
	df <- col-1							# 自由度
	o <- t(apply(dat, 1, rank))					# ケースごとに順位付けする
	R <- colSums(o)							# 条件ごとに順位の和をとる
	# chi <- 12*sum(R^2)/(row*col*(col+1))-3*row*(col+1)		# 検定統計量
	tie <- sum(apply(o, 1, function(x) {y <- table(x); sum(y^3-y)}))
	chi <- 12*sum((R-row*(col+1)/2)^2) / (row*col*(col+1)-tie/(col-1))
	p <- pchisq(chi, df, lower.tail=FALSE)				# P 値
	result1 <- structure(list(statistic=c("chi squared"=chi), 	# 検定結果のまとめ
		parameter=c(df=df), p.value=p, method=method,
		data.name=data.name), class="htest")
	R.m <- R/row							# 条件ごとの平均順位
	V <- sum((o-(col+1)/2)^2)					# 分散
	S <- combn(col, 2, function(ij) row^2*df*diff(R.m[ij])^2/(2*V))	# 対比較の検定統計量
	p <- pchisq(S, df, lower.tail=FALSE)				# P 値
	result2 <- cbind("chi sa."=S, "p-value"=p)
	rownames(result2) <- combn(col, 2, paste, collapse=":")
	return(structure(list(result1=result1, result2=result2), class="friedman"))
}
# print メソッド
print.friedman <- function(	obj,					# friedman が返すオブジェクト
				digits=4)				# 結果の表示桁数
{
	print(obj$result1, digits=digits)				# 全体として差があるかの検定結果
	cat("多重比較の結果\n\n")
	print(obj$result2, digits=digits)				# 多重比較の結果
}
# 分割表形式で与えられたデータに基づいて，グッドマン・クラスカルのガンマγを計算する
gamma2 <- function(f)				# 合計欄を含めない分割表
{
	R <- nrow(f)+2
	C <- ncol(f)+2
	g <- matrix(0, nr=R, nc=C)
	g[2:(R-1), 2:(C-1)] <- f
	cc <- dd <- 0
	for (i in 2:(R-1)) {
		for (j in 2:(C-1)) {
			cc <- cc+g[i,j]*sum(g[1:(i-1), 1:(j-1)], g[(i+1):R, (j+1):C])
			dd <- dd+g[i,j]*sum(g[1:(i-1), (j+1):C], g[(i+1):R, 1:(j-1)])
		}
	}
	return((cc-dd)/(cc+dd))
}
# 特定の相関係数行列を持つ多変量データの生成
gendat <- function(	nc,			# 標本サイズ
			r)			# 相関係数行列
{
	nv <- ncol(r)				# 変数の個数
	z <- matrix(rnorm(nv*nc), ncol=nv)	# 仮のデータ行列を作る。この時点では変数間の相関は近似的に0である。
	r2 <- cor(z)				# 相関係数行列
	res <- eigen(r2)			# 主成分分析を行い，
	coeff <-  solve(r2) %*% (sqrt(matrix(res$values, nv, nv, byrow=TRUE))*res$vectors)
	z <- t((t(z)-colMeans(z))/sqrt(apply(z, 2, var)*(nc-1)/nc)) %*% coeff	#主成分得点を求める。この時点で変数間の相関は完全に0である。
	return(z %*% chol(r))			# コレスキー分解の結果をもとに，指定された相関係数行列 を持つように主成分得点を変換する。
}
# 指定した平均値 mu と標準偏差 sigma を持つ n 個の正規乱数を発生させる
gendat1 <- function(n, mu=0, sigma=1)
{
	x <- rnorm(n)				# n 個の正規乱数を発生
	return((x-mean(x))/sd(x)*sigma+mu)	# 基準化する
}
# 特定の相関係数を持つ二変量データの生成
gendat2 <- function(	nc,						# サンプルサイズ
			r)						# 相関係数
{
	z <- matrix(rnorm(2*nc), ncol=2)				# 仮のデータ行列を作る。この時点では変数間の相関は近似的に0である。
	r2 <- cor(z)							# 変数間の相関係数を求める
	res <- eigen(r2)						# 主成分分析を行い，
	coeff <-  solve(r2) %*% t(sqrt(res$values)*t(res$vectors))	# 主成分得点係数を求める。
	z <- scale(z) %*% coeff						# 主成分得点を求める。この時点で変数間の相関は完全に0である。
	return(z %*% chol(matrix(c(1, r, r, 1), ncol=2)))		# コレスキー分解の結果をもとに，指定された相関係数行列 を持つように主成分得点を変換する。
}
# 一般化固有値問題を解く　Ax=λBx, A:実対称行列，B:実対称正値行列，λ:スカラー，x:列行列
geneig <- function(	a,				# 行列 A
			b)				# 行列 B
{
    a <- as.matrix(a)
    b <- as.matrix(b)
    if (nrow(a) == 1) {					# 1 行 1 列の場合
        res <- list(values=as.vector(a/b), vectors=as.matrix(1))
    }
    else {
        res <- eigen(b)
        g <- diag(1/sqrt(res$values))
        v <- res$vectors
        res <- eigen(g %*% t(v) %*% a %*% v %*% g)
        res$vectors <-v %*% g %*% res$vectors
    }
    return(res)
}
# 適合度の検定（exact test）
gft <- function(o,								# 度数
		p=rep(1/length(o), length(o)))					# 理論比
{

	printf <- function(fmt, ...)						# C 言語の printf をシミュレートする
	{
		cat(sprintf(fmt, ...))
	}

	gen_tab <- function(y)							# 度数表の発生
	{
		if (y == 1) {
			x2 <- sum((o-e)^2/e)
			if (x2 >= x0 || abs(x2-x0) <= EPSILON) {
				w2 <<- w2+exp(temp-sum(fact[o+1]))*prod(p^o)
			}
		}
		else {
			gen_tab(y-1)
			while (o[1]) {
				o[y] <<- o[y]+1
				o[1] <<- o[1]-1
				gen_tab(y-1)
			}
			o[1] <<- o[1]+o[y]
			o[y] <<- 0
		}
	}

	EPSILON <- 1e-10
	total <- sum(o)								# サンプルサイズ
	p <- p/sum(p)								# 理論比
	e <- p*total								# 期待値
	x0 <- sum((o-e)^2/e)							# カイ二乗値
	n <- length(o)								# カテゴリー数
	p_val <- pchisq(x0, n-1, lower.tail=FALSE)				# P 値
	printf("カイ二乗値は %g，自由度は %i，P値は %g\n", x0, n-1, p_val)	# 結果出力
	fact <- lfactorial(0:(total+1))						# 対数を取った階乗表
	temp <- fact[total+1]
	w2 <- 0
	o <- numeric(n)
	o[1] <- total
	gen_tab(n)								# 正確検定
	printf("正確なP値は %g\n", w2)
}
# 2 群のヒストグラム
hist2 <- function(	x1,							# 第一群のデータ
			x2,							# 第二群のデータ
			brks=NULL,						# 階級分割点
			...)							# barplot に引き渡す任意の引数
{
	if (is.null(brks)) {							# 階級分割点が与えられないときには，適切に設定
		brks <- hist(c(x1, x2), right=FALSE, plot=FALSE)$breaks
	}
	c1 <- hist(x1, breaks=brks, right=FALSE, plot=FALSE)$counts		# 度数1
	c2 <- hist(x2, breaks=brks, right=FALSE, plot=FALSE)$counts		# 度数2
	barplot(rbind(c1, c2), beside=TRUE, space=c(0, 0.2),			# 棒を並べて描く
		names.arg=brks[-length(c1)], axisnames=TRUE, axis.lty=1, ...)	# 横軸の目盛りラベル等
}
#####
#
# 独立 k 標本の平均値，標準偏差を求め，必要なら平均値・代表値の差の検定を行う
#
#####

indep.sample <- function(i,							# 分析対象の変数が入っているデータフレーム上の列番号または変数名ベクトル
			 j,							# 群を表す変数が入っているデータフレーム上の列番号または変数名ベクトル
			 df,							# データフレーム
			 latex=TRUE,						# LaTeX 形式で集計表を出力する（デフォルトは LaTeX 形式）
                         captions=NULL,						# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
                         labels=NULL,						# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
			 test=c("none", "parametric", "non-parametric"),	# デフォルト none では検定を行わない。検定を行うときはその種類を指定する
			 statistics=c("mean", "median"),			# （平均値，標準偏差）を求めるか（中央値，四分偏差）を求めるかを指定する
			 var.equal=FALSE,					# parametric の場合に等分散性を仮定するかどうかの引数
			 digits=3,						# 平均値，標準偏差を表示するときの小数点以下の桁数
			 output="",						# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
			 encoding=getOption("encoding"))			# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

# 下請け関数

	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

	SIQ <- function(x) return(diff(fivenum(x)[c(2,4)]))                     # 四分偏差を求める下請け関数

	indep.sample.sub <- function(ii, jj)
	{
		group <- colnames(df)[jj]					# 群を表す変数の名前
		df2 <- df[, c(ii, jj)]						# データフレームの列番号 ii, jj から 2 変数を取り出す
		df2 <- subset(df2, complete.cases(df2))				# 欠損値を持つケースを除く
		x <- df2[, 1]							# 分析対象変数
		g <- df2[, 2]							# 群変数
		lst <- list(g)							# by 関数を適用するために，群をリスト化する
		nt <- length(x)							# 全体のサンプルサイズ
		mt <- MEAN(x)							# 全体の平均値
		st <- SD(x)							# 全体の標準偏差
		n <- by(x, lst, length)						# 各群のサンプルサイズ
		m <- by(x, lst, MEAN)						# 各群の平均値（中央値）
		s <- by(x, lst, SD)						# 各群の標準偏差（四分偏差）
		nr <- length(table(lst))					# 群の数
		if (latex) {							# LaTeX 形式で集計結果を出力する
			cat("\n\\begin{table}[htbp]\n", file=output)		# \begin{table}[htbp]
			if (is.null(captions)) {
				cat(sprintf("\\caption{%s別の%sの集計}\n",	# \caption{○○別の□□の集計}
				    group, colnames(df2)[1]), file=output)
			}
			else {
				cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
			}
			if (!is.null(labels)) {
				cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
			}
			cat("\\centering\n", file=output)			# \centering
			cat("\\begin{tabular}{lccc} \\hline\n",file=output)	# \begin{tabular}{lccc} \hline
			cat(sprintf("& \\multicolumn{3}{c}{%s}\\\\ \\cline{2-4}\n", # \multicolumn{3}{c}{□□} \\ \cline{2-4}
			    colnames(df2)[1]), file=output)
			cat(group, "データ数", M.str, S.str, sep=" & ",		# ○○ & データ数 & 平均値 & 標準偏差
			    file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			for (l in 1:nr) {					# 各群について，
				cat(names(n)[l], n[l], sprintf(format, m[l]),	# 群の名前 & 数値 & 数値 & 数値
				    sprintf(format, s[l]), sep=" & ", file=output)
				cat("\\\\", file=output)			# \\
				if (l == nr) cat("\\hline\n", file=output)	# 最後の群なら \hline そうでなければ何もなし
				else cat("\n", file=output)
			}
			cat("全体", nt, sprintf(format, mt),			# 全体 & 数値 & 数値 & 数値
			    sprintf(format, st), sep=" & ", file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			cat("\\end{tabular}\n", file=output)			# \end{tabular}
		}
		else {								# tab で区切って出力する
			cat("\n表　", group, "別の", colnames(df2)[1],		# 表　○○別の□□の集計
			    "の集計", sep="", file=output)
			cat("\n", colnames(df2)[1], sep="\t", file=output,	# 　　□□
			    fill=TRUE)
			cat(group, "データ数", M.str, S.str, sep="\t",		# ○○　データ数　平均値　標準偏差
			    file=output, fill=TRUE)
			for (l in 1:nr) {					# 各群について，
				cat(names(n)[l], n[l], sprintf(format, m[l]),	# 群の名前　数値　数値　数値
				    sprintf(format, s[l]), sep="\t",
				    file=output, fill=TRUE)
			}
			cat("全体", nt, sprintf(format, mt),			# 全体　数値　数値　数値
			    sprintf(format, st), sep="\t",
			    file=output, fill=TRUE)
		}
		if (nr == 2) {							# 2 群の場合には，
			if (latex && test != "none") {				# LaTeX 形式で出力するときには，改行して次の行に出力準備
				cat("\\\\ \\noindent\n", file=output)
			}
			if (test == "parametric") {				# 平均値の差の検定のために t.test 関数を使う
				res <- t.test(x~g, var.equal=var.equal)		# t.test を呼ぶ
				cat(sprintf(if (latex) "$t$値 = %.3f, 自由度 = %.3f, $P$値 = %.3f\n"
					    else "t 値 = %.3f, 自由度 = %.3f, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
			else if (test == "non-parametric") {			# マン・ホイットニーの U 検定
				res <- wilcox.test(x~g)				# wilcox.test を呼ぶ
				cat(sprintf(if (latex) "$U$ = %.3f, $P$値 = %.3f\n"
					    else "U = %.3f, P 値 = %.3f\n",
					    res$statistic, res$p.value), file=output)
			}
		}
		else if (nr >= 3) {
			if (latex && test != "none") {				# LaTeX 形式で出力するときには，改行して次の行に出力準備
				cat("\\\\ \\noindent\n", file=output)
			}
			if (test == "parametric") {				# 一元配置分散分析
				res <- oneway.test(x ~ g, var.equal=var.equal)
				cat(sprintf(if (latex) "$F$値 = %.3f, 自由度 = (%i, %.3f), $P$値 = %.3f\n"
					    else "F 値 = %.3f, 自由度 = (%i, %.3f), P 値 = %.3f\n",
					    res$statistic, res$parameter[1], res$parameter[2], res$p.value), file=output)
			}
			else if (test == "non-parametric") {			# クラスカル・ウォリス検定
				res <- kruskal.test(x~g)	
				cat(sprintf(if (latex) "$\\chi_{kw}^2$ = %.3f, 自由度 = %i, $P$値 = %.3f\n"
					    else "カイ二乗値(kw) = %.3f, 自由度 = %i, P 値 = %.3f\n", 
					    res$statistic, res$parameter, res$p.value), file=output)
			} 
		}
		if (latex) {							# LaTeX 形式で集計結果を出力する場合は，
			cat("\\end{table}\n", file=output)			# \end{table}
		}
	}

# 関数本体
	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

	test <- match.arg(test)							# 引数が省略形で与えられたときに，正確な名前をとる
	statistics <- match.arg(statistics)					# 引数が省略形で与えられたときに，正確な名前をとる
	if (statistics == "mean") {
		MEAN <- mean							# 位置の母数を求める関数：平均値
		SD <- sd							# 散布度を求める関数：標準偏差
		M.str <- "平均値"
		S.str <- "標準偏差"
	}
	else {
		MEAN <- median							# 位置の母数を求める関数：中央値
		SD <-  SIQ							# 散布度を求める関数：四分偏差
		M.str <- "中央値"
		S.str <- "四分偏差"
	}
	format <- paste("%.", digits, "f", sep="")				# 小数点以下 3 桁で出力する書式
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}

	index <- 0
	for (jj in j) {								# j はベクトルまたはスカラーで，群を表す変数がある列番号
		for (ii in i) {							# i はベクトルまたはスカラーで，分析対象変数がある列番号
			if (ii != jj) {						# i, j の全ての組み合わせについて（ii と jj が違うときのみ），
				index <- index+1				# 集計のための下請け関数 indep.sample.sub を呼ぶ
				indep.sample.sub(ii, jj)
			}
		}
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
# 級内相関係数を計算する
intraclass.correlation <- function(dat)		# データ行列（n 行 2 列）
{
	dat <- subset(dat, complete.cases(dat))	# 欠損値を持つケースを除く
	nr <- nrow(dat)				# ケース数
	nc <- ncol(dat)				# 繰り返し数
	nt <- length(dat)			# データ総数					
	m <- rowMeans(dat)			# 行ごとの平均
	U <- apply(dat, 1, var)			# 行ごとの不偏分散
	gm <- mean(dat)				# 全体の平均
	Sw <- sum((nc-1)*U)			# 群内変動
	Mw <- Sw/(nt-nr)			# 群内分散
	Sb <- sum(nc*(m-gm)^2)			# 群間変動
	Mb <- Sb/(nr-1)				# 群間分散
	return((Mb-Mw)/(Mb+(nc-1)*Mw))		# 級内相関係数
}
# 相関係数のジャックナイフ推定
jack.knife <- function(x, y)					# 2 変数のベクトル
{
	OK <- complete.cases(x, y)				# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	n <- length(x)						# サンプルサイズ
	est <- cor(x, y)					# 標本相関係数
	cor.est <- sapply(1:n, function(i) cor(x[-i], y[-i]))	# i 組目のデータを除いたときの相関係数
	jk <- n*est-(n-1)*cor.est				# ジャックナイフ推定量
	plot(1:n, cor.est)					# 推定値の図を描いてみる
	list(sample.cor=est, jk.est=mean(jk))
}
# カッパ統計量
kappa.stat <- function(	o,						# 評定結果の二次元表
			w=FALSE,					# 重み行列（対角成分は何でも良い）
			conf.level=0.95)				# 信頼率
{
	data.name <- deparse(substitute(o))
	method <- "κ統計量"
	n <- sum(o)							# サンプルサイズ
	e <- outer(rowSums(o), colSums(o))/n				# 期待値
	if (is.matrix(w) == FALSE) {					# 重み付けしないときは，
		po <- sum(diag(o))/n
		qo <- 1-po
		pe <- sum(diag(e))/n
		qe <- 1-pe
		kappa <- 1-qo/qe
		sk <- sqrt(po*qo/(n*qe^2))				# κの信頼区間を計算するときに使う標準誤差
		sk0 <- sqrt(pe/(n*qe))					# H0：κ = 0 の検定に使うκの標準誤差
	}
	else {								# 重み付けするκを求めるときは，
		method <- paste(method, "（重み付け）", sep="")
		qo <- sum(w*o)/n
		qo2 <- sum(w*w*o)/n
		qe <- sum(w*e)/n
		qe2 <- sum(w*w*e)/n
		kappa <- 1-qo/qe
		sk <- sqrt((qo2-qo^2)/n/qe^2)				# κの信頼区間を計算するときに使う標準誤差
		sk0 <- sqrt((qe2-qe^2)/n/qe^2)				# H0：κ = 0 の検定に使うκの標準誤差
	}
	stopifnot(kappa >= 0)
	z <- kappa/sk0							# H0：κ = 0 の検定
	P <- pnorm(z, lower.tail=FALSE)*2				# P 値
	names(kappa) <- "kappa"
	names(z) <- "Z"
	cint <- kappa+c(-sk, sk)*qnorm(0.5+conf.level/2)
	attr(cint, "conf.level") <- conf.level
	return(structure(list(estimate=kappa, statistic=z, p.value=P,	# 結果は htest クラスのオブジェクト
		method=method, data.name=data.name, conf.int=cint, "sigma-kappa"=sk,
		"sigma-kappa-0"=sk0), class="htest"))
}
kendall2 <- function(x, y)		# 2変数間のケンドールの順位相関係数
{
	cor(cbind(x, y), use="complete.obs", method="kendall")[1, 2]
}

kendall <- function(dat)		# ケンドールの順位相関係数行列
{

	cor(dat, use="complete.obs", method="kendall")
}
# ケンドールの一致度係数
kendall.w <- function(x)				# データ行列
{
	method <- "ケンドールの一致度係数"
	data.name <- deparse(substitute(x))
	x <- subset(x, complete.cases(x))		# 欠損値を持つケースを除く	
	nv <- ncol(x)					# データ行列の列数（変数の個数）
	nc <- nrow(x)					# データ行列の行数（ケース数）
	o <- apply(x, 2, rank)				# 各列に rank 関数を適用（列ごとに順位を付ける）
	t <- apply(o, 2, table)				# 同順位を取る個数
	s1 <- rowSums(o)
	tie <- sapply(t, function(i) sum(i^3-i))	# 同順位の調整
	w <- 12*sum((s1-sum(s1)/nc)^2)/(nv^2*(nc^3-nc)	# ケンドールの W
		-nv*sum(tie))
	chi <- nv*(nc-1)*w				# 検定統計量（カイ二乗分布に従う）
	p <- pchisq(chi, nc-1, lower.tail=FALSE)	# P 値	
	return(structure(list(statistic=c("Kendall W"=w, "chi sq."=chi),
		parameter=c(df=nc-1), p.value=p, method=method,
		data.name=data.name), class="htest"))
}
# Kaplan-Meier 法による生命表
km.surv <- function(	time,						# 生存期間ベクトル
			event)						# エンドポイント（故障なら 1，そうでなければ 0）
{
	OK <- complete.cases(time, event)				# 欠損値を持つデータを除く
	time <- time[OK]
	event <- event[OK]
	fraction <- min(time[time > 0])/1000				# 打ち切りと故障の生存期間が同じときでも，ソートしたときに打ち切りの方が後になるような微小値を作る
	n <- length(time)						# 当初データ数
	truncate <- 1-event						# 打ち切りデータのとき 1 になるような変数を作る
	time <- time+truncate*fraction					# 打切りデータの生存期間に微小値を加える
	o <- order(time)						# 生存期間の順序づけ
	time <- time[o]							# 生存期間の並べ替え
	truncate <- truncate[o]						# 打ち切りかどうかのデータの並べ替え
	time <- time-truncate*fraction					# 生存期間データから微小値を取り除く
	i <- 1:n							# 整数ベクトル
	p <- ifelse(truncate, 1, (n-i)/(n-i+1))				# 生存確率は，打ち切りなら 1，そうでないなら (n-i)/(n-i+1)
	P <- cumprod(p)							# 累積生存率は累積
	se <- (1-truncate)/(n-i+1)/(n-i)				# 標準誤差の計算に使用する
	SE <- ifelse(truncate == 0, P*sqrt(cumsum(se)), NA)		# 標準誤差は故障時点でのみ計算される
	result <- data.frame(time, truncate, p, P, SE)			# 結果をデータフレームとして返す
	time <- c(0, time)						# 生存率グラフを描くために準備
	P <- c(1, P)							# 同じく
	plot(time, P, xlim=c(0, time[n+1]), ylim=c(0, 1), type="s")	# type="s" で階段状のグラフが描ける
	points(time, P)							# 打ち切りか故障があったポイントをマークする
	return(result)
}
kmo <- function(x)				# データ行列またはデータフレーム
{
	x <- subset(x, complete.cases(x))	# 欠損値を持つケースを除く
	r <- cor(x)                             # 相関係数行列
	r2 <- r^2                               # 相関係数行列の要素の二乗
	i <- solve(r)                           # 相関係数行列の逆行列
	d <- diag(i)                            # 対角成分
	p2 <- (-i/sqrt(outer(d, d)))^2          # 偏相関係数行列の要素の二乗
	diag(r2) <- diag(p2) <- 0               # 対角成分は計算には用いない
	KMO <- sum(r2)/(sum(r2)+sum(p2))
	MSA <- colSums(r2)/(colSums(r2)+colSums(p2))
	return(list(KMO=KMO, MSA=MSA))
}
# クラスカル・ウォリス検定（plus 多重比較）
kruskal.wallis <- function(	data,						# 測定値ベクトルまたはリスト
				group=NULL)					# 群変数ベクトル（data がリストの場合には無視される）
{
	method <- "クラスカル・ウォリス検定（plus 多重比較）"
	if (is.list(data)) {							# 第 1 引数がリストなら，
		data.name <- deparse(substitute(data))
		group <- factor(rep(1:length(data), sapply(data, length)))	# 群変数ベクトルを作る
		data <- unlist(data)						# リストをベクトルにする
	}
	else {
		data.name <- paste(deparse(substitute(data)), "~", deparse(substitute(group)))
	}
	OK <- complete.cases(data, group)					# 欠損値を持つケースを除く
	data <- data[OK]
	group <- group[OK]
	ni<- table(group)							# 各群のデータ数
	n <- sum(ni)								# 全データ数
	r <- rank(data)								# 順位付け
	Ri <- tapply(r, group, sum)						# 群ごとの和をとる
	S <- 12*sum(Ri^2/ni)/(n*(n+1))-3*(n+1)					# 検定統計量
	if (length(unique(data)) != n) {					# 同値があるなら
		tie <- table(data)
		S <- S/(1-sum(tie^3-tie)/(n^3-n))				# 同順位を考慮した検定統計量
	}
	a <- length(ni)								# 群の個数
	df <- a-1								# 自由度
	p <- pchisq(S, df, lower.tail=FALSE)					# P 値
	result1 <- structure(list(statistic=c("Kruskal-Wallis chi-squared"=S),
			parameter=c(df=df), p.value=p, method=method,
			data.name=data.name), class="htest")
	a.mean <- (n+1)/2							# 順位の平均値
	R.mean <- Ri/ni								# 順位和の平均値
	V <- sum((r-a.mean)^2)/(n-1)						# 分散
	S <- combn(a, 2, function(ij) diff(R.mean[ij])^2/sum(V/ni[ij]) )	# 一対比較
	p <- pchisq(S, df, lower.tail=FALSE)					# P 値
	result2 <- cbind("chi sq."=S, "p-value"=p)
	rownames(result2) <- combn(a, 2, paste, collapse=":")
	return(structure(list(result1=result1, result2=result2), class="kruskal.wallis"))
}
# print メソッド
print.kruskal.wallis <- function(	obj,					# kruskal.wallis が返すオブジェクト
					digits=4)				# 結果の表示桁数
{
	print(obj$result1, digits=digits)					# 全体として差があるかの検定結果
	cat("多重比較の結果\n\n")
	print(obj$result2, digits=digits)					# 多重比較の結果
}
# 二標本コルモゴロフ・スミルノフ検定を行う
ks2 <- function(	obs1,					# 第一群の度数分布
			obs2)					# 第二群の度数分布
{
	stopifnot(length(obs1) == length(obs2))			# 要素数は同じでなくてはならない
	name1 <- deparse(substitute(obs1))
	name2 <- deparse(substitute(obs2))
	data.name <- paste(name1, "and", name2)
	method <- "二標本コルモゴロフ・スミルノフ検定"
	n1 <- sum(obs1)						# 第一群のサンプルサイズ
	n2 <- sum(obs2)						# 第二群のサンプルサイズ
	cum1 <- cumsum(obs1)/n1					# 第一群の累積相対度数
	cum2 <- cumsum(obs2)/n2					# 第二群の累積相対度数
	D <- max(abs(cum1-cum2))				# 差の絶対値の最大値
	chi <- 4*D^2*n1*n2/(n1+n2)				# 検定統計量
	df <- 2
	p <- min(1, pchisq(chi, df, lower.tail=FALSE)*2)	# P 値
	stat <- c(D=D, "X-squared"=chi)
	names(df) <- "df"
	return(structure(list(statistic=stat, parameter=df, p.value=p,
		method=method, data.name=data.name,
		obs1=obs1, obs2=obs2, name1=name1, name2=name2),
		class=c("htest", "ks2")))			# 結果をまとめて返す
}
# plot メソッド（群別の度数分布図を描く）
plot.ks2 <- function(	obj,					# ks2 関数が返すオブジェクト
			method=c("barplot", "polygon"),		# barplot か度数多角形か
			name=c(obj$name1, obj$name2),		# 各群の名前
			col=c("grey30", "grey90"),		# barplot の色
			pch=c(19, 2),				# 度数多角形のマーク
			col2=c(1,4),				# 度数多角形の線の色
			lty=c(1,2),				# 度数多角形の線の種類
			xlab="", ylab="パーセント",		# 軸の名称
			x="topright", y=NULL,			# legend の位置
			...)					# barplot, matplot への引数
{
	d <- rbind(obj$obs1, obj$obs2)
	d <- d/rowSums(d)*100
	nc <- ncol(d)
	if (match.arg(method) == "barplot") {
		barplot(d, beside=TRUE, col=col, ylab=ylab, ...)
		legend(x, y, legend=name, fill=col)
		axis(1, at=1:nc*3-1, labels=1:nc, pos=0)
	}
	else {
		matplot(t(d), type="l", col=col2, lty=lty, xlab=xlab, ylab=ylab, xaxt="n", bty="n", ...)
		points(rep(1:nc, each=2), d, col=rep(col2, nc), pch=rep(pch, nc))
		legend(x, y, legend=name, col=col2, lty=lty, pch=pch)
		axis(1, at=1:nc, labels=1:nc, pos=0)
	}
}
# クラスカル・ウォリス検定（分割表データから）
kw3.test <- function(d)		# 分割表
{
	v <- rep(col(d), d)	# 測定値を再現
	g <- rep(row(d), d)	# 群変数を再現
	kruskal.test(v, g)	# R にある関数を呼び出す
}
# 集計表を与えて，クラスカル・ウォリス検定を行う
kw4.test <- function(tbl)				# 集計表（合計欄を含まない）
{

	method <- "集計表を与えてクラスカル・ウォリス検定"
	data.name <- deparse(substitute(tbl))
	nc <- ncol(tbl)					# 表頭に来る変数（順序変数）がとるカテゴリー数
	ni <- rowSums(tbl)				# 表側にくる変数が表す各群のデータ数
	n <- sum(tbl)					# 総合計
	tj <- colSums(tbl)				# 各カテゴリーに該当するデータ数
	rj <- c(0, cumsum(tj)[-nc])+(tj+1)/2		# 各カテゴリーに属するデータに与えられる順位
	Sx <- 12*sum((tbl%*%rj)^2/ni)/(n*(n+1))-3*(n+1)	# 検定統計量
	S0 <- Sx/(1-sum(tj^3-tj)/(n^3-n))		# 同順位を修正した検定統計量（カイ二乗分布に従う）
	df <- nrow(tbl)-1				# 自由度
	P <- pchisq(S0, df, lower.tail=FALSE)		# P 値
	return(structure(list(statistic=c("Kruskal-Wallis chi-squared"=S0), parameter=c(df=df), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# 共分散比を計算する
kyo.bunsan.hi <- function(dat)				# データ行列（合計点は含まない）
{
	dat <- subset(dat, complete.cases(dat))		# 欠損値を持つケースを除く
	total <- rowSums(dat)				# 合計点を計算する
	dat <- cbind(dat, total)			# 合計点の列を最右端列に加える
	nc <- ncol(dat)					# 列数
	SD <- apply(dat, 2, sd)				# 標準偏差を計算する
	result <- 100*(cor(dat)[nc,]*SD/SD[nc])[-nc]
	names(result) <- paste("Var", 1:(nc-1), sep="")
	return(result)
}
# 観察値と予測値の差の絶対値の p 乗の和を最小にするという方法による回帰直線のパラメータ推定
least.sum.abs <- function(	x,					# 独立変数ベクトル
				y,					# 従属変数ベクトル
				n=1,					# ブートストラップ法で信頼区間を求めるときの回数
				p=1,					# 1 のとき絶対値の和を最小にする。2 のときは普通の最小二乗法
				control=list())				# optim 関数への引数
{
	least.sum.abs0 <- function(x, y)				# 1 組のデータについて，切片と傾きの推定値を計算する
	{
		evaluate.error <- function(par)
		{
			return(sum(abs(y-par[2]*x-par[1])^p))
		}
		if (p < 1) {
			warning("p は 1 以上の値でなければなりません。p=1 として実行します。")
			p <- 1
		}
		ans <- optim(c(0, 0), evaluate.error, control=control)
		return(c(Intercept=ans$par[1], Slope=ans$par[2]))
	}
	Driver <- function(x, y)					# ブートストラップ法のためのドライバー
	{
		n <- length(x)
		suffix <- sample(n, n, replace=TRUE)			# リサンプリング
		return(least.sum.abs0(x[suffix], y[suffix]))		# リサンプリングしたデータについてパラメータを推定
	}
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	OK <- complete.cases(x, y)					# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	ans <- list(coefficients=least.sum.abs0(x, y),			# 引数に対してパラメータを推定する
		    names.xy=names.xy, x=x, y=y)
	if (n > 1) {
		ans2 <- replicate(n, Driver(x, y))			# ブートストラップを n 回実行
		ans <- append(ans, list(intercepts=ans2[1,], slopes=ans2[2,]))
	}
	class(ans) <- "least.sum.abs"					# print, plot メソッドのためにクラス名をつけておく
	return(ans)
}
# print メソッド
print.least.sum.abs <- function(obj,					# "least.sum.abs" オブジェクト
				digits=5,				# 表示桁数
				sig=0.95)				# 信頼度
{
	if (length(obj) == 4) {
		cat("Intercept:", round(obj$coefficients[1], digits),
		    "    Slope:", round(obj$coefficients[2], digits), "\n")
	}
	else {
		alpha <- (1-sig)/2
		LCL <- c(quantile(obj$intercepts, alpha), quantile(obj$slopes, alpha))
		UCL <- c(quantile(obj$intercepts, 1-alpha), quantile(obj$slopes, 1-alpha))
		ans <- data.frame(obj$coefficients, LCL, UCL)
		dimnames(ans) <- list(c("Intercept:", "Slope:"),
				      c("Estimate", paste(c(alpha, 1-alpha), "%", sep="")))
		print(round(ans, digits=digits))
	}
}
# plot メソッド
plot.least.sum.abs <- function(obj,					# "least.sum.abs" オブジェクト
			posx="topleft", posy=NULL,			# legend 関数のための位置引数
			xlab=obj$names.xy[1], ylab=obj$names.xy[2],	# 軸の名前
			hist=FALSE,					# ヒストグラムを描くとき TRUE にする
			...)						# その他の任意の plot 関数の引数
{
	if (hist && length(obj) == 6) {					# ブートストラップの結果を，hist=TRUE のときに，ヒストグラムで表示する
		layout(matrix(1:2, 2))
		hist(obj$intercepts, xlab="Intercept", main="", right=FALSE)
		hist(obj$slopes, xlab="Slope", main="", right=FALSE)
		layout(1)
	}
	else {								# 散布図と least.sum.abs 法の回帰直線と直線回帰式を表示する
		plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
		abline(obj$coefficients)
		abline(lm(obj$y~obj$x), lty=2, col=2)
		legend(posx, posy, legend=c("least sum abs", "linear regression"), lty=1:2, col=1:2)
	}
}
# 等分散性の検定
Levene.test <- function(x,					# データベクトル
			group,					# 群変数ベクトル
			method = c("mean", "median"))		# 検定統計量の計算方法
{
	data.name <- paste(deparse(substitute(x)), "~", deparse(substitute(group)))
	OK <- complete.cases(x, group)				# 欠損値を持つケースを除く
	x <- x[OK]
	fac <- as.factor(group[OK])
	fac <- fac[, drop=TRUE]
	n <- length(x)						# 全体のデータ個数
	n.i <- tapply(x, fac, length)				# 各群のデータ個数
	k <- length(n.i)					# 群の数
	method <- match.arg(method)				# 引数の補完
	x <- abs(x-tapply(x, fac, method)[fac])			# 測定値からそのデータが属する群の平均値（または中央値）を差し引く
	sw <- sum((n.i-1)*tapply(x, fac, var))			# 群内変動
	dfw <- n-k						# 群内変動の自由度
	dfb <- k-1						# 群間変動の自由度
	f <- ((var(x)*(n-1)-sw)/dfb)/(sw/dfw)			# 検定統計量
	P <- pf(f, dfb, dfw, lower.tail=FALSE)			# P 値
	method <- paste("等分散性の検定（", method, " で調整）", sep="")
	return(structure(list(statistic=c(F=f), parameter=c("df b"=dfb, "df w"=dfw), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# リッカート尺度を計算する
likert <- function(dat)							# 各カテゴリーへの回答数
{
	n <- length(dat)						# カテゴリー数
	resp <- dat/sum(dat)						# 各カテゴリーへの反応の割合（相対度数）
	cum <- cumsum(resp)						# 各カテゴリーへの反応の累積相対度数
	result <-(dnorm(qnorm(c(0, cum[-n])))-dnorm(qnorm(cum)))/resp
	names(result) <- paste("Cat", 1:n, sep="")			# 名前を付ける
	return(result)							# 結果を返す
}
# 重回帰分析における分散分析表
lm.anova <- function(obj) {
    df.reg <- obj$rank - 1
    df.res <- obj$df.residual
    f <- obj$fitted.values
    ss.reg <- sum((f - mean(f))^2)
    ss.res <- sum(obj$residuals^2)
    ms.reg <- ss.reg / df.reg
    ms.res <- ss.res / df.res
    F <- ms.reg / ms.res
    p.value <- pf(F, df.reg, df.res, lower.tail=FALSE)
    result <- data.frame(df=c(df.reg, df.res), SS=c(ss.reg, ss.res),
                         MS=c(ms.reg, ms.res), F=c(F, NA),
                         P=c(p.value, NA))
    colnames(result) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
    rownames(result) <- c("Regression", "Residuals")
    class(result) <- c("anova", "data.frame")
    return(result)
}
# 対数正規確率紙に累積相対度数をプロットする
lnpp <- function(x)					# データベクトル
{
	log.axis <- function(z)				# 対数目盛りの軸を描く
	{
		z <- floor(log10(z))			# 対数にしたときの整数部
		log.min <- min(z)			# 最小値
		z2 <- 1:10*10^log.min			# 一番左側に描かれる目盛り数値のセット
		n <- max(z)-log.min			# 10 倍しながら順次，右の位置に目盛りを描く
		z2 <- rep(z2, n+1)*10^rep(0:n, each=10)	# 対数目盛り位置の数値
		log.z2 <- log10(z2)			# 目盛りを描く位置
		axis(1, at=log.z2, labels=z2)		# log.z2 の位置に，z2 という数値を描く
		abline(v=log.z2, col="gray")		# 垂直格子線を入れる
	}
	n <- length(x)					# データの個数
	log.x <- log10(sort(x))				# データをプロットするときの横座標
	y <- ((1:n)-0.5)/n					# 累積確率
	probs <- c(0.01, 0.1, 1, 5, 10, 20, 30, 40, 50,	# 縦軸の目盛り
		   60, 70, 80, 90, 95, 99, 99.9, 99.99)/100
	plot(log.x[c(1,n)], qnorm(probs[c(1,17)]),	# 枠組み
		type="n", xaxt="n", yaxt="n",
		xlab="Observed Value", ylab="Cumulative Percent",
		main="Log Normal Probability Paper")
	log.axis(x)					# 横軸（対数目盛）を描く
	axis(2, qnorm(probs), probs*100)		# 縦軸（正規確立目盛）を描く
	abline(h=qnorm(probs), col="grey")		# 水平格子線を描く
	points(log.x, qnorm(y))				# データ点を描く
}
# 片対数軸または両対数軸でプロット（散布図）を描く
log.plot <- function(	x,				# 横軸に取る変数
			y,				# 縦軸に取る変数
			log.of.x=FALSE,			# 横軸も対数軸にするときに TRUE にする
			x.label="",			# 横軸の名前
			y.label="",			# 縦軸の名前
			title="",			# 図のタイトル
			color="gray")			# 格子線の色
{
	log.axis <- function(z, which)			# 対数軸を描く関数
	{
		z <- floor(log10(z))			# 対数にしたときの整数部
		log.min <- min(z)			# 最小値
		z2 <- 1:10*10^log.min			# 値の範囲をカバーするように
		n <- max(z)-log.min			# 10 倍しながら順次，右の位置に目盛りを描く
		z2 <- rep(z2, n+1)*10^rep(0:n, each=10)	# 対数目盛り位置の数値
		log.z2 <- log10(z2)			# 目盛りを描く位置
		axis(which, at=log.z2, labels=z2)	# log.z2 の位置に，z2 という数値を描く
		if (which == 1) {
			abline(v=log.z2, col=color)	# 垂直格子線を描く
		}
		else {
			abline(h=log.z2, col=color)	# 水平格子線を描く
		}
	}

	OK <- complete.cases(x, y)			# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	log.y <- log10(y)				# 縦軸に取る変数の常用対数を取る
	if (log.of.x == FALSE) {			# 縦軸だけが対数軸の場合
		plot(x, log.y, type="n", yaxt="n", xlab=x.label, ylab=y.label, main=title)
		points(x, log.y)			# データ点を描く
		abline(v=x, col=color)			# 垂直格子線を描く
	}
	else {						# 縦軸，横軸ともに対数軸の場合
	  	log.x <- log10(x)			# 横軸に取る変数の常用対数を取る
	 	plot(log.x, log.y, type="n", xaxt="n", yaxt="n", xlab=x.label, ylab=y.label, main=title)
	 	points(log.x, log.y)			# データ点を描く
		log.axis(x, 1)				# 横軸を対数軸として描く
	}
	log.axis(y, 2)					# 縦軸を対数軸として描く
}
# 多重ロジスティックモデル(Walker-Duncan 法)
logistic.regression <- function(data)							# データ行列（最右端の列が従属変数）
{

	printf <- function(fmt, ...)							# 書式付き print 関数
	{
		cat(sprintf(fmt, ...))
	}

	diff <- function(x, m, n, coeff)						# 数値微分
	{
		mp1 <- m+1
		temp <- coeff[mp1]+coeff[1:m]%*%t(x[, 1:m, drop=FALSE])
		p0 <- 1/(1+exp(-temp))
		p1 <- 1-p0
		pp <- ifelse(x[, mp1] == 1, p1, -p0)
		diff1 <- numeric(mp1)
		diff1[mp1] <- sum(pp)
		diff1[1:m] <- diff1[1:m]+colSums(x[,1:m, drop=FALSE]*pp)
		temp <- -x[, 1:m, drop=FALSE]*p0*p1
		diff2 <- matrix(0, nrow=mp1, ncol=mp1)
		diff2[1:m, 1:m] <- t(x[,1:m, drop=FALSE])%*%as.matrix(temp)
		diff2[lower.tri(diff2, diag=FALSE)] <- 0
		diff2[mp1, mp1] <- -sum(p0*p1)
		diff2[1:m, mp1] <- diff2[1:m, mp1]+colSums(temp[,1:m, drop=FALSE])
		diff2 <- diff2+t(diff2)
		diag(diff2) <- diag(diff2)/2
		return(list(diff1=diff1, diff2=diff2))
	}

	llh <- function(x, m, coeff)							# 対数尤度
	{
		temp <- colSums(t(x[,1:m])*coeff[1:m])+coeff[m+1]
		-sum(log(1+exp(ifelse(x[,m+1] == 1, -temp, temp))))
	}

	newton.logist <- function(x, m, n, sds)						# ニュートン法によるあてはめ
	{
		mp1 <- m+1
		coeff0 <- numeric(mp1)
		coeff <- rep(1e-14, mp1)
		for (itr in 1:500) {
			temp <- diff(x, m, n, coeff)
			coeff0 <- solve(temp$diff2, temp$diff1)
			converge <- all(abs(coeff0/coeff) < 1e-10)
			coeff <- coeff-coeff0
			if (converge) {
				break
			}
		}
		stopifnot(converge == TRUE)
		se <- sqrt(-diag(solve(temp$diff2)))
		printf("\n対数尤度 = %.14g\n", llh(x, m, coeff))
		printf("\nパラメータ推定値\n\n")
		printf("         %14s %14s %12s %8s %14s\n", "偏回帰係数", "標準偏差", "t 値", "P 値", "標準化偏回帰係数")
		t <- abs(coeff[1:m]/se[1:m])
		p <- pt(t, n-mp1, lower.tail=FALSE)*2
		for (i in 1:m) {
			printf("%8s %14.7g %14.7g %12.7g %8.5f %14.7g\n", vnames[i], coeff[i], se[i], t[i], p[i], coeff[i]*sds[i])
		}
		t <- abs(coeff[mp1]/se[mp1])
		p <- pt(t, n-mp1, lower.tail=FALSE)*2
		printf("%8s %14.7g %14.7g %12.7g %8.5f\n%60s%i\n\n", "　定数項", coeff[mp1], se[mp1], t, p, "t の自由度 = ", as.integer(n-mp1))
		return(coeff)
	}

	fitness <- function(data, m, n, coeff)						# 当てはまり具合の確認
	{
		lambda <- coeff[1:m]%*%t(data[,1:m])+coeff[m+1]				# λ
		pred <- 1/(1+exp(-lambda))						# リスクの予測値
		y <- data[,m+1]								# エンドポイント
		div <- round(seq(0, n-1, by=n/10),0)[-1]				# 対象をほぼ同数に10区分するときのサンプル数
		xs <- sort(lambda)							# λをソートして
		div2 <- (xs[div]+xs[div+1])/2						# 分割点を求める
		g <- findInterval(lambda, div2)+1					# λがどの区分にはいるか群分け
		from <- c(min(lambda), xs[div+1])					# λの区間 [from, to]
		to <- c(xs[div], max(lambda))
		mid <- (from+to)/2							# 区間の級中心
		pred <- tapply(pred, g, sum)						# エンドポイントを持つものの期待値
		obs <- tapply(y, g, sum)						# エンドポイントを持つものの観察値
		cnt <- tapply(y, g, length)						# サンプル数
		valid <- as.integer(names(pred))					# 区間幅が正のもの
		from <- from[valid]							# 調整
		to <- to[valid]								# 調整
		table <- data.frame("以上"=from, "以下"=to,
				"期待値"=pred, "リスク"=pred/cnt,
				"観察値"=obs, "故障率"=obs/cnt, "サンプル数"=cnt)
		print(table)
		printf("%60s%s\n", "", "左の2列は，各区間のλの値（最小値と最大値）")
		plot(c(min(from), max(to)), c(0, max(ifelse(pred>obs, pred, obs)/cnt)),	# リスク，故障率，理論曲線の描画
			type="n", xlab="Lambda", ylab="Risk")
		for (i in 1:10) {
			lines(c(from[i], to[i]), rep(pred[i]/cnt[i], 2), lty=3)
			lines(c(from[i], to[i]), rep(obs[i]/cnt[i], 2))
		}
		x <- seq(min(lambda), max(lambda), length.out=100)
		lines(x, 1/(1+exp(-x)))
	}

# logistic 関数本体

	data <- subset(data, complete.cases(data))					# 欠損値を除く
	n <- nrow(data)									# サンプルサイズ
	mp1 <- ncol(data)								# 列数
	m <- mp1-1									# 独立変数の個数
	vnames <- colnames(data)							# 変数名を取り出す
	if (is.null(vnames)) vnames <- paste("変数", 1:mp1, sep="")			# 変数名がなかったら定義する
	num <- table(data[, mp1])							# サンプルの内訳

	printf("***** 多重ロジスティック回帰\n\n")
	printf("サンプルサイズ　　　　%5i\n", n)
	printf("　　生存（打ち切り）　%5i\n", as.integer(num[1]))
	printf("　　死亡（故障）　　　%5i\n", as.integer(num[2]))
	if (num[1] == 0 || num[2] == 0 || num[1]+num[2] == 2 || n <= mp1) {
		stop("有効ケース数が 1 以下です\n")
	}

	means <- colMeans(data)								# 各変数の平均値
	sds <- apply(data, 2, sd)							# 各変数の標準偏差
	printf("\n          %15s %15s\n", "平均値", "標準偏差")	
	for (i in 1:m) {
		printf("%8s  %15.7g %15.7g\n", vnames[i], means[i], sds[i])
	}

	coeff <- newton.logist(data, m, n, sds)						# パラメータ推定
	fitness(data, m, n, coeff)							# 当てはまり具合の確認
	invisible(coeff)								# 係数のみを返す（自動的に表示はされない）
}
# ログランク検定を行う
logrank <- function(	group,					# 群を識別するベクトル（1, 2 のいずれか）
			event,					# 死亡なら 1，生存なら 0 の値をとるベクトル
			time,					# 生存期間ベクトル
			method=c("SAS", "Tominaga"))		# 一般的な SAS などの方法か，富永の方法か
{
	method <- match.arg(method)
	data.name <- sprintf("time: %s, event: %s, group: %s",
		deparse(substitute(time)), deparse(substitute(event)), deparse(substitute(group)))
	OK <- complete.cases(group, event, time)		# 欠損値を持つケースを除く
	group <- group[OK]
	event <- event[OK]
	time <- time[OK]
	len <- length(group)
	stopifnot(length(event) == len, length(time) == len)
	tg <- table(c(time, rep(NA, 4)),			# 後ろの4項目はダミー
		    c(group, 1, 1, 2, 2)*10+c(event, 1, 0, 1, 0))
	k <- nrow(tg)
	nia <- table(group)[1]
	nib <- len-nia
	na <- c(nia, (rep(nia, k)-cumsum(tg[,1]+tg[,2]))[-k])
	nb <- c(nib, (rep(nib, k)-cumsum(tg[,3]+tg[,4]))[-k])
	da <- tg[,2]
	db <- tg[,4]
	dt <- da+db
	nt <- na+nb
	d <- dt/nt
	O <- c(sum(da), sum(db))
	ea <- na*d
	eb <- nb*d
	E <- c(sum(ea), sum(eb))

	result <- data.frame(da, db, dt, na, nb, nt, d, ea, eb)

	if (method == "Tominaga") {				# 富永による検定方式
		method <- "ログランク検定（富永）"
		chi <- sum((O-E)^2/E)
	}
	else {							# SAS 等と同じ検定方式
		method <- "ログランク検定（一般的）"
		v <- sum(dt*(nt-dt)/(nt-1)*na/nt*(1-na/nt), na.rm=TRUE)
		chi <- (sum(da)-sum(na*d))^2/v
	}
	P <- pchisq(chi, 1, lower.tail=FALSE)
	return(structure(list(statistic=c("X-squared"=chi), parameter=c(df=1), p.value=P,
	method=method, data.name=data.name, result=result), class="htest"))
}
# 多倍長計算により，n 番目のフィボナッチ数を返す
longFibonacci <- function(n)
{
	"%add%" <- function(ans, b)                             # 足し算の演算子　ans %add% b を行い結果を返す
	{                                                       # ans, b は多倍長整数
		if (length(ans) != length(b)) {
			ans <- c(ans, 0)
		}
		ans <- ans+b                                    # 各桁の足し算を行う
		if (ans[length(ans)] >= 10000000000) {
			ans <- c(ans, 0)
		}
		for (i in 1:length(ans)) {                      # 各桁について下の桁から，
			if (ans[i] >= 10000000000) {            # 繰り上がり処理を行う
				ans[i] <- ans[i]-10000000000
				ans[i+1] <- ans[i+1]+1
			}
		}
		return(ans)                                     # 結果を返す
	}
	if (n <= 2) {
		c <- 1
	}
	else {
		a <- b <- c <- numeric(1)
		a[1] <- b[1] <- 1
		for (i in 3:n) {
			c <- a %add% b
			a <- b
			b <- c
		}
	}
	class(c) <- "longFibonacci"
	return(c)
}
print.longFibonacci <- function(x)				# プリント・メソッド
{
	top.zero <- "          "
	for (i in length(x):1) {
		if (x[i] != 0 || top.zero == FALSE) {
			out <- paste(top.zero, as.character(x[i]), sep="")
			len <- nchar(out)
			cat(sprintf("%10s ", substring(out, len-9, len)))
			top.zero <- "0000000000"
		}
	}
	cat("\n")
}
# ライアンの方法とチューキーの方法による平均値の対比較
m.multi.comp <- function(	n,				# 標本サイズベクトル
				me,				# 平均値ベクトル
				s,				# 標準偏差ベクトル
				alpha=0.05,			# 有意水準
				method=c("ryan", "tukey"))	# 方法
{
	printf <- function(fmt, ...)
	{
		cat(sprintf(fmt, ...))
	}

	check <- function(s, b)					# 検定しようとしている二群が，それまでに有意でないとされた二群に挟まれているか
	{
		if (ns.n > 1) {
			for (i in 1:ns.n) {
				if (ns.s[i] <= s && s <= ns.b[i] && ns.s[i] <= b && b <= ns.b[i]) {
					return(FALSE) 		# 検定するまでもなく有意でないとする
				}
			}
		}
		return(TRUE)					# 検定しなくてはならない
	}

	k <- length(n)						# 群の数
	stopifnot(k == length(me), k == length(s), n > 0, floor(n) == n, s > 0)
	method <- match.arg(method)				# 引数の補完
	o <- order(me)						# 平均値の大きさの順位
	sn <- n[o]						# 並べ替えた標本サイズ
	sm <- me[o]						# 並べ替えた平均値
	ss <- s[o]						# 並べ替えた標準偏差
	nt <- sum(sn)						# 全体の標本サイズ
	mt <- sum(sn*sm)/nt					# 全体の平均値
	dfw <- nt-k						# 群内平方和の自由度
	vw <- sum((sn-1)*ss^2)/dfw				# 群内分散
	num.significant <- ns.n <- 0
	ns.s <- ns.b <- numeric(k*(k-1)/2)			# 有意でない群の記録用
	for (m in k:2) {					# 検定対象の選定
		for (small in 1:(k-m+1)) {
			big <- small+m-1
			if (check(small, big)) {
				t0 <- (sm[big]-sm[small])/sqrt(vw*(1/sn[big]+1/sn[small]))	# 検定統計量
				if (method == "ryan") { 					# Ryan の方法
					P <- pt(t0, dfw, lower.tail=FALSE)*2			# 有意確率
					nominal.alpha <- 2*alpha/(k*(m-1))			# 名義的有意水準
					result <- P <= nominal.alpha				# 検定結果
				}
				else { # Tukey の方法
					t0 <- t0*sqrt(2)
					q <- (qtukey(0.05, k, dfw, lower.tail=FALSE) + qtukey(0.05, m, dfw, lower.tail=FALSE))/2
					WSD <- q*sqrt(vw*(1/sn[big]+1/sn[small])/2)		# Wholly Significant Difference
#					P <- ptukey(t0, m, dfw, lower.tail=FALSE)		# 有意確率　以下で置き換え
					P <- uniroot(function(p)((qtukey(p, k, dfw, lower.tail=FALSE)+
								  qtukey(p, m, dfw, lower.tail=FALSE))/2 -t0), c(0,1),tol=1e-7)$root
					result <- sm[big]-sm[small] >= WSD			# 検定結果
				}
				if (result) { 			# 有意であるとき
					num.significant <- 1
					o.small <- o[small]	# 群番号は，入力時（実験時）の順序番号を表示する
					o.big <- o[big]
					printf("mean[%2i]=%7.5f vs. mean[%2i]=%7.5f : diff.= %7.5f, ",
						o.small, me[o.small], o.big, me[o.big], me[o.big]-me[o.small])
					if (method == "ryan") {
						printf("t=%7.5f : P=%7.5f, alpha'=%7.5f\n", t0, P, nominal.alpha)
					}
					else {
						printf("WSD=%7.5f : t=%7.5f : P=%7.5f\n", WSD, t0, P)
					}
				}
				else {				# 有意でないとき
					ns.n <- ns.n+1
					ns.s[ns.n] <- small
					ns.b[ns.n] <- big
				}
			}
		}
	}
	if (num.significant == 0) {				# 有意差のある群は一つもなかった
		print("Not significant at all.")
	}
}
# マハラノビスの距離による基準群への帰属確率
Mahalanobis <- function(dat,			# 基準群のデータ行列
			x)			# 所属確率を計算するデータ行列
{
	dat <- subset(dat, complete.cases(dat))	# 欠損値を持つケースを除く
	n <- nrow(dat)				# ケース数
	p <- ncol(dat)				# 変数の個数
	ss <- var(dat)*(n-1)/n			# 分散・共分散行列
	inv.ss <- solve(ss)			# 分散共分散行列の逆行列
	m <- colMeans(dat)			# 各変数の平均値
	dif <- t(t(x)-m)			# 平均値からの偏差
	d2 <- apply(dif, 1,
		function(z) z %*% inv.ss %*% z) # マハラノビスの平方距離
	P <- pchisq(d2, p, lower.tail=FALSE)	# 所属確率
	return(data.frame(d2=d2, P=P))
}
# アイテムデータをカテゴリーデータに変換する
make.dummy <- function(dat)
{
	ncat <- ncol(dat)
	dat[, 1:ncat]  <- lapply(dat, function(x) {
		if (is.factor(x)) {
			return(as.integer(x))
		}
		else {
			return(x)
		}
	})
	mx <- sapply(dat, max)
	start <- c(0, cumsum(mx)[1:(ncat-1)])
	nobe <- sum(mx)
	t(apply(dat, 1, function(obs) 1:nobe %in% (start+obs))) + 0
}
# 塗り分け地図を描く
map <- function(code.list,						# 描画する都道府県コード
		density=NULL,						# ハッチングの 1 インチあたりの線密度
		color=NULL)						# 塗り分けに使用する色
{
	map0 <- function(data, dens, color)				# 描画関数
	{
		continue <- apply(data, 1, any)				# 経度と緯度が共に 0 が描画の区切り
		plot(lon, lat, type = "n", axes=FALSE,			# 大枠を決める
			xlab="", ylab="", bty="n", asp=1)
		start <- 1
		k <- 0
		for (i in 2:nrow(data)) {
			if (continue[i] == FALSE) {			# 区切れ目
				k <- k+1				# 何番目の描画か
				if (i-start == 4) {			# 沖縄を描くときの区画線
					lines(data[start:(i-1),])
				}
				else {
					polygon(data[start:(i-1),], density=dens[k], col=color[k], border="black")
				}
				start <- i+1
			}
		}
	}

# 関数本体

	for (i in seq(along=code.list)) {				# 指定した全ての都道府県について
		if (code.list[i] %in% c(15, 28, 47)) {			# 新潟県，兵庫県，沖縄県は描画パーツが 2 つ
			code.list <- c(code.list, -code.list[i])	# コードリストの追加（目印とするために負の数で）
			density <- c(density, density[i])		# 線密度と
			color <- c(color, color[i])			# 色も追加する
		}
	}
	code.list[code.list == -15] <- 48				# 追加のコードリストに直す
	code.list[code.list == -28] <- 49
	code.list[code.list == -47] <- 50

	lon <- lat <- NULL
	for (i in code.list) {						# 指定した全ての都道府県について
		fn <- sprintf("jpn/%02i", i)				# データファイル名を得る
		gwm <- matrix(scan(fn, quiet=TRUE), ncol=2, byrow=TRUE)	# データを読む
		lon <- c(lon, gwm[,1], 0)				# 経度を蓄積
		lat <- c(lat, gwm[,2], 0)				# 緯度を蓄積
	}
	mlon <- min(lon[lon != 0])					# 経度の最小値
	mlat <- max(lat[lat != 0])					# 緯度の最大値
	lon <- ifelse(lon == 0, 0, lon-mlon+1)				# 経度の調整
	lat <- ifelse(lat == 0, 0, mlat-lat+1)				# 緯度の調整
	map0(cbind(as.integer(lon), as.integer(lat)), density, color)	# 描画関数を呼ぶ
}
# 度数分布表から中央値を求める
median2 <- function(	f,					# 度数ベクトル
			b,					# 最初の階級の下限値
			w)					# 階級幅
{
	cf <- cumsum(f)						# 累積度数
	n <- sum(f)						# 標本サイズ
	position <- length(cf[cf < n/2])			# 中央値の存在する階級
	b+w*position+w*(n/2-cf[position])/f[position+1]
}
# 指定した平均値ベクトル mu と標準偏差ベクトル sigma を持つ母集団から，比率ベクトル prob で抽出される n 個の混合正規乱数を発生させる
mix1 <- function(	n,					# データ数
			mu,					# 平均値ベクトル
			sigma,					# 標準偏差ベクトル
			prob=rep(1/length(mu),length(mu)))	# 抽出割合
{
	k <- length(mu)
	if (k == length(sigma) && k == length(prob)) {
		suffix <- sample(k, n, replace=TRUE, prob=prob)
		x <- mapply(function(mean, sd) rnorm(n, mean, sd), mu, sigma)
		return(list(d=x[cbind(1:n,suffix)], which=suffix))
	}
}
# 指定した平均値ベクトル mu と標準偏差ベクトル sigma を持つ母集団から，比率ベクトル prob で抽出される n 個の混合正規乱数を発生させる
mix2 <- function(	n,					# データ組数
			mu,					# 平均値行列（グループ×2変数）
			sigma,					# 標準偏差行列（グループ×2変数）
			r,					# 相関係数ベクトル（グループ）
			prob=rep(1/length(mu),length(mu)))	# 抽出割合（グループ）
{
	library(MASS)
	k <- nrow(mu)
	if (k == nrow(sigma) && k == length(r) && k == length(prob) &&
	    ncol(mu) == 2 && ncol(sigma) == 2) {
		suffix <- sample(k, n, replace=TRUE, prob=prob)
		x <- mapply(function(mean1, mean2, sd1, sd2, r) {
			cat(mean1, mean2, sd1, sd2, r, "\n")
			mvrnorm(n, mu=c(mean1, mean2),
			        Sigma=matrix(c(sd1^2, r*sd1*sd2, r*sd1*sd2, sd2^2), 2),
			        empirical=TRUE)},
			mu[,1], mu[,2], sigma[,1], sigma[,2], r)
		dim(x) <- c(n, 2, k)
		return(list(d=t(mapply(function(i, k) x[i,,k], 1:n, suffix)), which=suffix))
	}
}
# 牧厚志ら「経済・経営のための統計学」有斐閣アルマの
# 第8章「マーケティング ブランド選択行動の分析」濱岡豊
# 本文に示されている R プログラムを一般化して使いやすくした
mlm <- function(	d,			# 本文中で $d_{ij}$ と記されているもの
			x=NULL,		# 本文中で $x_{ijm}$ と記されているもの（変数の順序（列）に注意）
			z=NULL,		# 本文中で $z_{i}$ と記されているもの
			constants=TRUE)	# 本文中で「ブランド定数」と記されているものを加えるか否か
{						# x, z, constant は 1 つ以上選択可。途中の引数が省略される場合には引数名付きで指定すること
	getV <- function(par)				# 効用 $V_{ij}$ の計算
	{
		V <- matrix(0, n, N)
		if (include.x) {				# $x_{ijm}$ について
			beta <- par[1:p]
			V <- V+t(apply(x3, 3, function(z) z%*%beta))
		}
		if (include.z) {				# $z_{i}$ について
			beta <- cbind(matrix(par[include.x*p+1:(ncol.z*(n-1))], ncol.z, n-1), rep(0, ncol.z))
			V <- V+t(z%*%beta)
		}
		if (constants) {				# ブランド定数について
			V <- V+c(par[(n.parameters-n+2):n.parameters], 0)
		}
		return(V)
	}

	mlogit <- function(par)				# 対数尤度を求める関数
	{
		V <- getV(par)
		return(sum(V[cbind(d, 1:N)]-log(colSums(exp(V)))))
	}

	include.x <- !is.null(x)				# x が指定されると TRUE
	include.z <- !is.null(z)				# z が指定されると TRUE
	stopifnot(include.x || include.z || constants)	# どれか 1 つは指定されていないといけない
	if (is.data.frame(d)) {				# データフレームなら行列にしておく
		d <- data.matrix(d)
	}
	n <- length(unique(d))				# 選択対象の種類
	N <- length(d)					# サンプルサイズ
	if (include.x && is.data.frame(x)) {		# データフレームなら行列にしておく
		x <- data.matrix(x)
		p <- ncol(x)/n				# 選択対象の属性の種類
		if (p != floor( p )) {
			stop("属性の総項目数が選択対象の種類の整数倍になっていない")
		}
		x3 <- array(unlist(x), dim=c(N, p, n))	# 選択対象別に計算するために 3 次元配列にする
	}
	if (include.z && is.data.frame(z)) {		# データフレームなら行列にしておく
		z <- data.matrix(z)
	}
	n.parameters <- 0					# 推定すべきパラメータ数
	if (include.x) {
		n.parameters <- n.parameters+p		# + 選択対象の属性の種類
	}
	if (include.z) {
		ncol.z <- ncol(z)				# + z として使う変数の個数*(選択対象の種類-1)
		n.parameters <- n.parameters+ncol.z*(n-1)
	}
	if (constants) {
		n.parameters <- n.parameters+n-1		# + 選択対象の種類-1
	}
	par <- numeric(n.parameters)			# パラメータベクトル（初期値は 0 でよい）
	res <- optim(par, mlogit, method="BFGS", hessian=TRUE, control=list(fnscale=-1))
	LL <- res$value					# 対数尤度
	beta.estimated <- res$par				# 推定されたパラメータ
	AIC <- -2*(LL-n.parameters)				# AIC
	ses <- sqrt(diag(solve(-res$hessian)))		# 標準誤差
	t <- abs(beta.estimated)/ses			# t 値
	P <- pt(t, N-n.parameters, lower.tail=FALSE)*2	# P 値
	df <- data.frame(beta.estimated, ses, t=round(t, 3), P=round(P, 3))
	if(include.x) {					# 名前の編集
		rownames(df)[1:p] <- c(sub("[0-9]*$", "", colnames(x)[1:p]))
	}
	if (include.z) {
		for (i in 1:ncol.z) {
			rownames(df)[include.x*p++(i-1)*(n-1)+1:(n-1)] <- paste(colnames(z)[i], 1:(n-1), sep="-")
		}
	}
	if (constants) {
		rownames(df)[(n.parameters-n+2):n.parameters] <- paste("Constant", 1:(n-1), sep="-")
	}
	V <- t(getV(beta.estimated))			# 効用
	v <- exp(V)
	P <- v/rowSums(v)					# 選択確率
	predict <- apply(P, 1, which.max)			# 選択予測（P が一番大きいもの）
	res <- list(df=df, LL=LL, AIC=AIC, V=V, P=P, predict=predict)
	class(res) <- "mlm"
	return(res)
}

summary.mlm <- function(object, digits=5)			# sumamry メソッド
{
	print(object$df, digits=digits)
	cat("LL  =", object$LL, "\nAIC =", object$AIC, "\n")
}

print.mlm <- function(object, digits=5)			# print メソッド
{
	print.default(object, digits=digits)
}

predict.mlm <- function(object)				# predict メソッド
{
	object$predict
}
# 重回帰分析
mreg <- function(	dat,					# データ行列
			func.name=c("solve", "ginv"))		# 逆行列を計算する関数の選択
{
	dat <- subset(dat, complete.cases(dat))			# 欠損値を持つケースを除く
	n <- nrow(dat)						# ケース数
	nc <- ncol(dat)						# 列数
	nv <- nc-1						# 独立変数の個数（最後の一つは従属変数）
	if (is.null(colnames(dat))) {				# 変数名が付いていないときには仮の名前を付ける
		colnames(dat) <- paste("Var", 1:nc, sep="")
	}
	r <- cor(dat)						# 相関係数行列
	m <- colMeans(dat)[-nc]					# 独立変数の平均値ベクトル
	if (match.arg(func.name) == "solve") {
		inverse <- solve
		betas <- inverse(r[-nc, -nc], r[,nc][-nc])	# 標準化偏回帰係数
	}
	else {
		library(MASS)
		inverse <- ginv
		betas <- inverse(r[-nc, -nc]) %*% r[,nc][-nc]	# 標準化偏回帰係数
	}
	variance <- var(dat)*(n-1)				# 変動共変動行列
	prop <- diag(variance)					# 対角成分
	prop <- (prop/prop[nc])[-nc]				# 偏回帰係数に変換するための係数（独立変数の変動/従属変数の変動）
	b <- betas/sqrt(prop)					# 偏回帰係数
	Sr <- variance[,nc][-nc]%*%b				# 回帰による変動
	St <- variance[nc, nc]					# 全変動
	Se <- St-Sr						# 誤差変動
	SS <- c(Sr, Se, St)					# 平方和（変動）ベクトル
	dfr <- nv						# 回帰による変動の自由度
	dfe <- n-nv-1						# 誤差変動の自由度
	dft <- n-1						# 全変動の自由度
	df <- c(dfr, dfe, dft)					# 自由度ベクトル
	Ms <- SS/df						# 平均平方ベクトル
	f <- Ms[1]/Ms[2]					# F 値
	fvalue <- c(f, NA, NA)
	p <- c(pf(f, dfr, dfe, lower.tail=FALSE), NA, NA) 	# P 値
	b0 <- mean(dat[,nc])-sum(b*m)				# 定数項
	b <- c(b, b0)						# 偏回帰係数ベクトル
	inv <- inverse((n-1)*cov(dat)[-nc, -nc])
	SEb <- c(sapply(1:nv, function(i) sqrt(inv[i, i]*Ms[2])), sqrt((1/n+m%*%inv%*%m)*Ms[2]))
	tval <- b/SEb						# 偏回帰係数の有意性検定
	pval <- pt(abs(tval), n-nv-1, lower.tail=FALSE)*2	# P 値
	tolerance <- 1/diag(inverse(cor(dat)[-nc, -nc]))	# トレランス
	result <- cbind(b, SEb, tval, pval, c(betas, NA), c(tolerance, NA))
	rownames(result) <- c(colnames(dat)[1:nv], "定数項")
	colnames(result) <- c("偏回帰係数", "標準誤差", "t 値", "P 値", "標準化偏回帰係数", "トレランス")
	R2 <- 1-Se/St						# 重相関係数の二乗
	R <- sqrt(R2)						# 重相関係数
	R2s <- 1-Ms[2]/Ms[3]					# 自由度調整済み重相関係数の二乗
	loglik <- -0.5*n*(log(2*pi)+1-log(n)+log(Se))	# 対数尤度
	AIC <- 2*nc+2-2*loglik					# AIC
	Rs <- c("重相関係数"=R, "重相関係数の二乗"=R2, "自由度調整済重相関係数の二乗"=R2s,
		"対数尤度"=loglik, "AIC"=AIC)
	anova <- cbind(SS, df, Ms, fvalue, p)			# 分散分析表
	rownames(anova) <- c("回帰", "残差", "全体")
	colnames(anova) <- c("平方和", "自由度", "平均平方", "F 値", "P 値")
	return(structure(list(result=result, anova=anova, Rs=Rs), class="mreg"))
}
# print メソッド
print.mreg <- function(	obj,					# mreg が返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	print(obj$result, digits=digits, na.print="")
	cat("\n回帰の分散分析表\n\n")
	print(obj$anova, digits=digits, na.print="")
	cat("\n")
	sapply(1:length(obj$Rs), function(i) cat(names(obj$Rs)[i], "=", round(obj$Rs[i], digits), "\n"))
}
# 重回帰分析（mreg 関数を 3 カ所書き換えた）
mreg2 <- function(dat)						# データ行列
{
	dat <- subset(dat, complete.cases(dat))			# 欠損値を持つケースを除く
	n <- nrow(dat)						# ケース数
	nc <- ncol(dat)						# 列数
	nv <- nc-1						# 独立変数の個数（最後の一つは従属変数）
	if (is.null(colnames(dat))) {				# 変数名が付いていないときには仮の名前を付ける
		colnames(dat) <- paste("Var", 1:nc, sep="")
	}
	r <- cor(dat)						# 相関係数行列
	m <- colMeans(dat)[-nc]					# 独立変数の平均値ベクトル
#	betas <- solve(r[-nc, -nc], r[,nc][-nc])		# 標準化偏回帰係数
	betas <- ginv(r[-nc, -nc]) %*% r[,nc][-nc]		# 上の 1 行をこれに置き換え
	variance <- var(dat)*(n-1)				# 変動共変動行列
	prop <- diag(variance)					# 対角成分
	prop <- (prop/prop[nc])[-nc]				# 偏回帰係数に変換するための係数（独立変数の変動/従属変数の変動）
	b <- betas/sqrt(prop)					# 偏回帰係数
	Sr <- variance[,nc][-nc]%*%b				# 回帰による変動
	St <- variance[nc, nc]					# 全変動
	Se <- St-Sr						# 誤差変動
	SS <- c(Sr, Se, St)					# 平方和（変動）ベクトル
	dfr <- nv						# 回帰による変動の自由度
	dfe <- n-nv-1						# 誤差変動の自由度
	dft <- n-1						# 全変動の自由度
	df <- c(dfr, dfe, dft)					# 自由度ベクトル
	Ms <- SS/df						# 平均平方ベクトル
	f <- Ms[1]/Ms[2]					# F 値
	fvalue <- c(f, NA, NA)
	p <- c(pf(f, dfr, dfe, lower.tail=FALSE), NA, NA) 	# P 値
	b0 <- mean(dat[,nc])-sum(b*m)				# 定数項
	b <- c(b, b0)						# 偏回帰係数ベクトル
#	inv <- solve((n-1)*cov(dat)[-nc, -nc])	
	inv <- ginv((n-1)*cov(dat)[-nc, -nc])			# 上の 1 行をこれに置き換え
	SEb <- c(sapply(1:nv, function(i) sqrt(inv[i, i]*Ms[2])), sqrt((1/n+m%*%inv%*%m)*Ms[2]))
	tval <- b/SEb						# 偏回帰係数の有意性検定
	pval <- pt(abs(tval), n-nv-1, lower.tail=FALSE)*2	# P 値
#	tolerance <- 1/diag(solve(cor(dat)[-nc, -nc]))		# トレランス
	tolerance <- 1/diag(ginv(cor(dat)[-nc, -nc]))		# 上の 1 行をこれに置き換え
	result <- data.frame(b, SEb, tval, pval, c(betas, NA), c(tolerance, NA))
	rownames(result) <- c(colnames(dat)[1:nv], "定数項")
	colnames(result) <- c("偏回帰係数", "標準誤差", "t 値", "P 値", "標準化偏回帰係数", "トレランス")
	R2 <- 1-Se/St						# 重相関係数の二乗
	R <- sqrt(R2)						# 重相関係数
	R2s <- 1-Ms[2]/Ms[3]					# 自由度調整済み重相関係数の二乗
	loglik <- 0.5*(sum(-n*(log(2*pi)+1-log(n)+log(Se))))	# 対数尤度
	AIC <- 2*nc+2-2*loglik					# AIC
	Rs <- c("重相関係数"=R, "重相関係数の二乗"=R2, "自由度調整済重相関係数の二乗"=R2s,
		"対数尤度"=loglik, "AIC"=AIC)
	anova <- data.frame(SS, df, Ms, fvalue, p)		# 分散分析表
	rownames(anova) <- c("回帰", "残差", "全体")
	colnames(anova) <- c("平方和", "自由度", "平均平方", "F 値", "P 値")
	return(list(result=result, anova=anova, Rs=Rs))
}
#####
#
# マルチアンサー項目と別の 1 つの変数についてクロス集計表を作成する
#
#####

ma <- function(i,								# 表側に来る変数が入っているデータフレーム上の列番号または変数名ベクトル
	       j,								# 表側に来る変数が入っているデータフレーム上の列番号または変数名ベクトル
										# i, j いずれかがマルチアンサー項目を表すのでベクトルになる
	       df,								# データフレーム
	       latex=TRUE,							# LaTeX 形式で度数分布表を出力する（デフォルトは LaTeX 形式）
               caption=NULL,							# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
               label=NULL,							# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
	       output="",							# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
	       encoding=getOption("encoding"))					# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

# 下請け関数
	is.error <- function(k)							# k 列目の変数が factor であり，かつ，二値変数でなければならない
	{
		if (!is.factor(df[,k])) {					# factor でない
			warning(sprintf("%s を factor にしてください", colnames(df)[k]))
			return(TRUE)
		}
		else if (nlevels(df[,k]) != 2) {				# 二値変数でない
			warning(sprintf("%s が二値データではありません", colnames(df)[k]))
			return(TRUE)
		}
		return(FALSE)
	}

	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

# 関数本体
	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}
	if (length(i) == 1 && length(j) > 1) {					# i の方が普通の変数（1 つ），j の方がマルチアンサー変数（複数）のとき
		row <- TRUE							# 表頭にマルチアンサー項目，row = TRUE として識別
	}
	else if (length(i) > 1 && length(j) == 1) {				# i の方がマルチアンサー変数（複数），j  の方が普通の変数（1 つ）のとき
		row <- FALSE							# 表側にマルチアンサー項目，row = FALSE として識別
		temp <- i							# row == TRUE のときと同じように処理するために，i, j の内容を入れ替える
		i <- j
		j <- temp
	}
	else {									# 普通の変数が複数個指定されるような場合は想定していない
		stop("このプログラムでは，i, j いずれかが要素 1 のスカラー，他方が要素2以上のベクトルであることを仮定しています。使用法に誤りがあります")
	}
	if (any(sapply(j, is.error))) stop("データの準備に問題があります")	# マルチアンサー変数のチェック

	df.i <- as.factor(df[,i])						# 普通の変数の方も factor にされている方が望ましいが，そうではないときのために
	ans <- sapply(j, function(k) table(df[,i], df[,k]))[1:nlevels(df.i),]	# マルチアンサー集計表本体を作る
	ans <- cbind(ans, table(df.i))						# 合計列を作る（普通の変数の方の度数分布）
	ans <- rbind(ans, colSums(ans))						# 合計行を作る
	rownames(ans) <- c(levels(df.i), "合計")				# 行の名前
	colnames(ans) <- c(colnames(df[,j]), "該当数")				# 列の名前
	nc <- ncol(ans)								# 列数
	pct <- ans/ans[,nc]*100							# 行方向パーセント
	if (!row) {								# row == FALSE の場合に，
		ans <- t(ans)							# ans を転置
		pct <- t(pct)							# pct を転置
		nc <- ncol(ans)							# 転置後の列数を再計算
	}
	nr <- nrow(ans)								# 最終時点の集計表の行数
	if (latex) {								# LaTeX 形式で集計結果を出力する
		cat("\\begin{table}[htbp]\n", file=output)			# \begin{table}[htbp] 
		if (is.null(caption)) {
			cat("\\caption{マルチアンサー項目の集計}\n", file=output) # \caption{マルチアンサー項目の集計}
		}
		else {
			cat(sprintf("\\caption{%s}\n", caption), file=output)	# 指定した表題
		}
		if (!is.null(label)) {
			cat(sprintf("\\label{%s}\n", label), file=output)	# 指定したラベル
		}
		cat("\\centering\n", file=output)									# \centering
		cat("\\begin{tabular}{l", rep("c", nc), "} \\hline\n", sep="", file=output)				# \begin{tabular}{cc…c} \hline
		if (row) {												# 表頭にマルチアンサー項目
			cat(colnames(df)[i], colnames(ans), sep=" & ", file=output)					# 変数名 & マルチアンサー変数i & …
			cat("\\\\ \\hline\n", file=output)								# \\ \hline
			for (i in 1:nr) {										# 各行について，
				cat(rownames(ans)[i], ans[i,], sep=" & ", file=output)					# 行名 & 集計数i & …
				cat("\\\\\n", file=output)								# \\
				cat("\\%", sprintf("{\\small \\textit{%.1f}}", pct[i,-nc]), sep=" & ", file=output)	# % & パーセントi & …
				cat("\\\\", file=output)								# \\
				if (i >= nr-1) cat("\\hline\n", file=output) else cat("\n", file=output)		# 最終行の前なら \\
			}
		}
		else {													# 表側にマルチアンサー項目
			cat(sprintf("& \\multicolumn{%i}{c}{%s}\\\\ \\cline{2-%i}\n", nc-1, colnames(df)[i], nc),	# マルチアンサーではない方の変数名
			    file=output)
			cat("", colnames(ans), sep=" & ", file=output)							# 列名（変数のカテゴリー名）
			cat("\\\\ \\hline\n", file=output)								# \\ \hline
			for (i in 1:nr) {										# 各行（マルチアンサー項目）ごとに，
				cat(rownames(ans)[i], ans[i,], sep=" & ", file=output)					# 行名 & 集計数i & …
				cat("\\\\\n", file=output)								# \\
				if (i < nr) {										# 最終行でないときは，
					cat("\\%", sprintf("{\\small \\textit{%.1f}}", pct[i,]), sep=" & ", file=output)# % & パーセントi & …
					cat("\\\\", file=output)							# \\
				}
				if (i >= nr-1) cat("\\hline\n", file=output) else cat("\n", file=output)		# 最終行の前なら \\
			}
		}
		cat("\\end{tabular}\n", file=output)									# \end{tabular}
		cat("\\end{table}\n", file=output)									# \end{table}
	}
	else {									# tab で区切って出力する
		cat("表　マルチアンサー項目の集計\n", file=output)							# マルチアンサー項目の集計
		if (row) {
			cat("\n", file=output, fill=TRUE)
			cat(colnames(df)[i], colnames(ans), sep="\t", file=output, fill=TRUE)				# 変数名　マルチアンサー変数i …
			for (i in 1:nr) {										# 各行について，
				cat(rownames(ans)[i], ans[i,], sep="\t", file=output, fill=TRUE)			# 行名　　集計数i　　　…
				cat("%", sprintf("%.1f", pct[i,-nc]), sep="\t", file=output, fill=TRUE)			# % 　　　パーセントi　…
			}
		}
		else {													# 表側にマルチアンサー項目を置く場合
			cat("\n", colnames(df)[i], sep="\t", file=output, fill=TRUE)					# マルチアンサーではない方の変数名
			cat("", colnames(ans), sep="\t", file=output, fill=TRUE)					# 列名（変数のカテゴリー名）
			for (i in 1:nr) {										# 各行（マルチアンサー項目）ごとに，
				cat(rownames(ans)[i], ans[i,], sep="\t", file=output, fill=TRUE)			# 行名　集計数i　…
				if (i < nr) {										# 最終行でないときは，
					cat("%", sprintf("%.1f", pct[i,]), sep="\t", file=output, fill=TRUE)		# % 　　パーセントi　…
				}
			}
		}
	}
	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
# 多倍長計算を行う
"%add%" <- function(ans, b)				# 足し算の演算子　ans %add% b を行い結果を返す
{							# ans, b は "multibyte" クラスの多倍長整数
	ans <- ans+b					# 各桁の足し算を行う
	for (i in length(ans):1) {			# 各桁について下の桁から，
		if (ans[i] >= 10000000000) {		# 繰り上がり処理を行う
			ans[i] <- ans[i]-10000000000
			ans[i-1] <- ans[i-1]+1
		}
	}
	return(ans)					# 結果を返す
}
#
"%sub%" <- function(ans, b)				# 引き算の演算子　ans %sub% b を行い結果を返す
{							# ans, b は "multibyte" クラスの多倍長整数
	ans <- ans-b					# 各桁の引き算を行う
	for (i in length(ans):1) {			# 各桁について下の桁から，
		if (ans[i] < 0) {			# 繰り下がり処理を行う
			ans[i] <- ans[i]+10000000000
			ans[i-1] <- ans[i-1]-1
		}
	}
	return(ans)					# 結果を返す
}
#
"%div%" <- function(ans, n)				# 割り算の演算子　ans %div% n を行い結果を返す
{							# 注：n は "multibyte" クラスではなく普通の整数値
	r <- 0						# 剰余
	for (i in 1:length(ans)) {			# 各桁について上の桁から，
		x <- ans[i]+r*10000000000		# より上の位での剰余を考慮した，被除算数
		ans[i] <- x%/%n				# 割り算を行い結果を格納
		r <- x-n*ans[i]				# 今回の剰余
	}
	return(ans)					# 結果を返す
}
#
calc.pi <- function(len)				# πの計算例　小数点以下 len 桁まで求める
{
	len <- len %/% 10+3				# "multibyte" クラスの多倍長整数の必要桁数
	pi <- a <- b <- numeric(len)			# 多倍長変数の領域確保　pi <- a <- b <- 0
	a[1] <- 16*5					# a[1] が最上位桁（a は小数点以下の数値を格納）
	n <- 1						# 除数
	repeat {					# 十分な精度を持つまで繰り返し計算
		a <- a %div% 25				# a <- a / 25
		b <- a %div% n				# b <- a / n
		pi <- pi %add% b			# pi <- pi + b
		n <- n+2				# n <- n + 2
		a <- a %div% 25				# a <- a / 25
		b <- a %div% n				# b <- a / n
		pi <- pi %sub% b			# pi <- pi - b
		n <- n+2				# n <- n + 2
		if (sum(b) == 0) break 			# pi の値が変化しなくなったらループ終了
	}
	a <- b <- numeric(len)				# 多倍長変数の領域確保　a <- b <- 0
	a[1] <- 4*239					# a <- 4*239
	n <- 1						# 除数
	repeat {					# 十分な精度を持つまで繰り返し計算
		a <- a %div% 57121			# a <- a / 57121
		b <- a %div% n				# b <- a / n
		pi <- pi %sub% b			# pi <- pi - b
		n <- n+2				# n <- +2
		a <- a %div% 57121			# a <- a / 57121
		b <- a %div% n				# b <- a / n
		pi <- pi %add% b			# pi <- pi + b
		n <- n+2				# n <- n + 2
		if (sum(b) == 0) break			# pi の値が変化しなくなったらループ終了
	}
	class(pi) <- "multibyte"			# "multibyte" クラスにする（プリント・メソッドを使うため）
	return(pi)					# 結果を返す
}
#
calc.e <- function(len)					# e の計算例　小数点以下 len 桁まで求める
{
	len <- len %/% 10+3				# "multibyte" クラスの多倍長整数の必要桁数
	e <- t <- numeric(len)				# 多倍長変数の領域確保　e <- t <- 0
	e[1] <- t[1] <- 1				# e <- t <- 1
	i <- 0						# 除数
	repeat {					# 十分な精度を持つまで繰り返し計算
		i <- i+1				# i = 1, 2, ...
		t <- t %div% i				# t = 1 / i!
		if (sum(t) == 0) break			# t が 0 になるまで
		e <- e %add% t				# e <- e + t
	}
	class(e) <- "multibyte"				# "multibyte" クラスにする（プリント・メソッドを使うため）
	return(e)					# 結果を返す
}
#
print.multibyte <- function(ans)			# プリント・メソッド
{
	if (ans[1] == 3) {				# πの計算結果
		cat("π = 3.\n")
	}
	else {						# e の計算結果
		cat("e = 2.\n")
	}
	for (i in 2:(length(ans)-2)) {			# 小数点以下の答えを出力する（ちょっと冗長だが）
		out <- paste("0000000000", as.character(ans[i]), sep="")
		len <- nchar(out)
		cat(sprintf(" %10s", substring(out, len-9, len)))
		if ((i-1) %% 5 == 0) cat("\n")		# 1行に50桁ずつ出力する
	}
	cat("\n")
}
# 多項分布の確率を求める
multinomial <- function(x,							# 観察値ベクトル
			p)							# 母比率ベクトル
{
	stopifnot(length(x) == length(p), sum(p) == 1, floor(x) == x, x >= 0)
	exp(lgamma(sum(x)+1)+sum(x*log(p))-sum(sapply(x+1, lgamma)))		# 対数で計算して後で逆対数を求める
}

# 重相関係数を計算する
multiple.cor <- function(x)			# データ行列
{
	x <- subset(x, complete.cases(x))	# 欠損値を持つケースを除く
	r <- sqrt(1-1/diag(solve(cor(x))))	# その変数とその変数以外の変数の重相関係数を，変数ごとに計算する
	var.names <- colnames(x)
	names(r) <- if (is.null(var.names)) paste("Var", 1:ncol(x)) else var.names
	return(r)
}
# カイ二乗分布を用いる独立性の検定
my.chisq.test <- function(x)				# 分割表
{
	if (is.matrix(x)) {
		nr <- nrow(x)
		nc <- ncol(x)
		if (nr < 2 || nc < 2) {
			stop("行数・列数は 2 以上でないといけません")
		}
	}
	else {
		stop("分割表を入力してください")
	}
	data.name <- deparse(substitute(x))
	method <- "カイ二乗分布を用いる独立性の検定（残差分析）"
	rt <- rowSums(x)
	ct <- colSums(x)
	n <- sum(x)
	expectation <- outer(rt, ct)/n
	if (any(expectation < 1)) {
		warning("expectation less than 1")
	}
	else if (sum(expectation <= 5)/(nr*nc) > 0.2) {
		warning("more than 20% cells have expectation less than 5")
	}
	e <- (x-expectation)/sqrt(expectation)
	d <- e/sqrt(outer(1-rt/n, 1-ct/n))
	chi2 <- sum(e^2)
	df <- (nr-1)*(nc-1)
	P <- pchisq(chi2, df, lower.tail=FALSE)
	names(chi2) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=chi2, parameter=df, p.value=P,
		method=method, data.name=data.name, observed=x, expected=expectation,
		standardized.residuals=e, adjusted.residuals=d),
		class=c("htest", "my.chisq.test")))
}
# summary メソッド
summary.my.chisq.test <- function(	obj,		# my.chisq.test が返すオブジェクト
					digits=5)	# 出力桁数
{
	cat("調整された残差\n")
	print(obj$adjusted.residuals, digits=digits)
	cat("\nP 値\n")
	print(pnorm(abs(obj$adjusted.residuals), lower.tail=FALSE)*2, digits=digits)
}
# バートレット検定（分散の均一性の検定）
my.bartlett.test <- function(	n,					# 各群のデータ個数のベクトル
				u)					# 各群の不偏分散のベクトル
{
	stopifnot(	length(n) == length(u),				# 入力データのチェック
			n > 1,
			u > 0,
			floor(n) == n)
	data.name <- sprintf("%s for sample sizes and %s for variances",
		deparse(substitute(n)), deparse(substitute(u)))
	method <- "バートレット検定（分散の均一性の検定）"
	ng <- length(n)							# 群の個数
	temp1 <- n-1
	temp2 <- sum(temp1)
	chi.sq0 <- temp2*log(sum(temp1*u)/temp2) - sum(temp1*log(u))	# 検定統計量
	co <- 1+(sum(1/temp1)-1/temp2)/(3*ng-3)				# 修正項
	chi.sq <- chi.sq0/co						# 検定統計量（カイ二乗分布に従う）
	df <- ng-1
	p <- pchisq(chi.sq, df, lower.tail=FALSE)			# P 値
	return(structure(list(statistic=c("X-squared"=chi.sq),
		parameter=c(df=df), p.value=p, method=method,
		data.name=data.name), class="htest"))
}
# 正準相関分析
my.cancor <- function(	x,						# データ行列
			gr1,						# 第一群の変数がある列位置のベクトル
			gr2)						# 第二群の変数がある列位置のベクトル
{
	geneig2 <- function(a, b, k, sd)				# 一般化固有値問題を解く関数
	{
	    a <- as.matrix(a)
	    b <- as.matrix(b)
	    if (nrow(a) == 1) {
	        res <- list(values=a/b, vectors=as.matrix(1))
	    }
	    else {
	        res <- eigen(b)
	        g <- diag(1/sqrt(res$values))
	        v <- res$vectors
	        res <- eigen(g %*% t(v) %*% a %*% v %*% g)
	        res$vectors <-v %*% g %*% res$vectors
	    }
	    std.vectors <- res$vectors[,1:k]
	    unstd.vectors <- std.vectors/sd
	    list(values=res$values[1:k], std.vectors=std.vectors, unstd.vectors=unstd.vectors)
	}

	k <- min(length(gr1), length(gr2))				# 第一変数群と第二変数群の個数の少ない方の個数
	x <- subset(x, complete.cases(x))				# 欠損値を持つケースを除く
	r <- cor(x)							# 相関係数行列
	S11 <- r[gr1, gr1, drop=FALSE]					# 第一変数群の相関係数
	S22 <- r[gr2, gr2, drop=FALSE]					# 第二変数群の相関係数
	S12 <- r[gr1, gr2, drop=FALSE]					# 第一変数群と第二変数群の相関係数
	x1 <- as.matrix(x[, gr1, drop=FALSE])				# 第一変数群のデータ行列
	x2 <- as.matrix(x[, gr2, drop=FALSE])				# 第二変数群のデータ行列
	sd1 <- apply(x1, 2, sd)						# 第一変数群の標準偏差
	sd2 <- apply(x2, 2, sd)						# 第二変数群の標準偏差
	res1 <- geneig2(S12 %*% solve(S22) %*% t(S12), S11, k, sd1)	# 第一変数群に対する解
	res2 <- geneig2(t(S12) %*% solve(S11) %*% S12, S22, k, sd2)	# 第二変数群に対する解
	score1 <- scale(x1 %*% res1[[3]])				# 第一変数群に対する正準得点
	score2 <- scale(x2 %*% res2[[3]])				# 第二変数群に対する正準得点
	list(canonical.correlation.coefficients=sqrt(res1[[1]]),
		standardized.coefficients=list(group1=res1[[2]], group2=res2[[2]]),
		coefficients=list(group1=res1[[3]], group2=res2[[3]]),
		canonical.scores=list(group1=score1, group2=score2))
}
# 単相関係数，偏相関係数，重相関係数を計算する
my.cor <- function(x)				# データ行列
{
	x <- subset(x, complete.cases(x))	# 欠損値を持つケースを除く
	r <- cor(x)				# 相関係数行列
	i <- solve(r)				# 逆行列
	d <- diag(i)				# 対角要素
	p <- -i/sqrt(outer(d, d))		# 偏相関係数行列
	r[lower.tri(r)] <- p[lower.tri(p)]	# 単相関係数行列の下三角行列を偏相関係数行列にする
	diag(r) <- sqrt(1-1/d)			# 対角要素を重相関係数に置き換える
	var.names <- colnames(x)
	rownames(r) <- colnames(r) <- if (is.null(var.names)) paste("Var", 1:ncol(x)) else var.names
	return(r)
}
# 2×2 分割表のフィッシャーの正確確率検定
my.fisher <- function(x)					# 2×2 分割表
{
	odds.ratio <- function(a, b, c, d)			# オッズ比を求める関数
	{
		if (a*b*c*d == 0) {				# セルのどれかが 0 の場合
			a <- a+0.5				# それぞれに 0.5 を加える
			b <- b+0.5
			c <- c+0.5
			d <- d+0.5
		}
		res <- a*d/(b*c)				# オッズ比
		return(max(res, 1/res))				# どちらか大きい方を返す
	}

	stats <- function(i) {
		e <- (a <- i[1]) + (b <- i[2])
		f <- (c <- i[3]) + (d <- i[4])
		n <- (g <- a+c) + (h <- b+d)
		return(c(n*(a*d-b*c)^2/(e*f*g*h),		# カイ二乗値
			 odds.ratio(a, b, c, d),		# オッズ比
			 exp(lchoose(e, a)+lchoose(f, c)-lchoose(n, g))))	# 生起確率
	}

	ct <- colSums(x)					# 列和
	rt <- rowSums(x)					# 行和
	n <- sum(x)						# 総合計
	mx <- min(rt[1], ct[1])					# a が取り得る最大値
	mi <- max(0, rt[1]+ct[1]-n)				# a が取り得る最小値
	A <- mi:mx						# a のベクトル
	B <- rt[1]-A						# b のベクトル
	C <- ct[1]-A						# c のベクトル
	D <- ct[2]-B						# d のベクトル
	Cell <- cbind(A, B, C, D)				# 行方向の 4 つの数値が一つの分割表
	result<-apply(Cell, 1, stats)				# 行方向の 4 つの数値に，stats 関数を適用
	Chi.sq <- result[1,]					# カイ二乗値ベクトル
	OddsRatio <- result[2,]					# オッズ比ベクトル
	Probability <- result[3,]				# 生起確率ベクトル
	Cum.Prob1 <- cumsum(Probability)			# 生起確率の累積和
	Cum.Prob2 <- rev(cumsum(rev(Probability)))		# 逆方向からの生起確率の累積和
	Pearson <- Chi.sq >= stats(c(x))[1]-1e-15		# カイ二乗値から判定した P 値に算入すべきセル
	Fisher <- Probability <= stats(c(x))[3]+1e-15		# フィッシャーの定義による P 値に算入すべきセル
	OR <- OddsRatio >= stats(c(x))[2]-1e-15			# オッズ比から判定した P 値に算入すべきセル
	p.Pearson <- sum(Probability[Pearson])			# カイ二乗値から判定した P 値
	p.Fisher <- sum(Probability[Fisher])			# フィッシャーの定義による P 値
	p.OR <- sum(Probability[OR])				# オッズ比から判定した P 値
	return(list(result=cbind(Cell, Chi.sq, Probability, OddsRatio, Pearson, Fisher, OR, Cum.Prob1, Cum.Prob2), p.Pearson=p.Pearson, p.Fisher=p.Fisher, p.OR=p.OR))
}
# 中央値付近に同値がある場合に比例配分によりより妥当な中央値を求める
my.median <- function (	x,		# データベクトルまたは級限界のベクトル
			y=NULL,		# x がデータベクトルの場合は NULL，
					# x が級限界のベクトルの場合は度数ベクトル
			accuracy=0)	# 測定精度。省略時は普通のメディアン
{
	median.sub <- function(x)
	{
		n <- length(x)
		half <- (n + 1)/2
		if (n%%2 == 1) {
			sort(x, partial = half)[half]
		}
		else {
			sum(sort(x, partial = c(half, half + 1))[c(half, half + 1)])/2
		}
	}

	if (is.null(y)) {
		x <- x[!is.na(x)]
		median <- median.sub(x)
		if ((ntie <- length(x[x == median])) != 1 && accuracy > 0) {
			x <- c(x[x != median], (median-(ntie+1)*accuracy/2/ntie)+1:ntie/ntie)
			median <- median.sub(x)
		}
		return(median)
	}
	else {
		stopifnot(length(x)-1 == length(y))
		k <- length(y)
		csum <- cumsum(y)
		n <- csum[k]
		for (i in 1:k) {
			if (csum[i] >= n/2) break
		}
		return(x[i]-accu/2+(n/2-csum[i-1])/y[i]*(x[i+1]-x[i]))
	}
}
# 二次データに基づき，一元配置分散分析を行う
my.oneway.ANOVA <- function(	n,							# サンプルサイズ
				m,							# 平均値
				u,							# 不偏分散
				var.equal=FALSE)					# 等分散を仮定するとき TRUE
{
	data.name <- sprintf("%s for sample sizes, %s for means and %s for variances",
				deparse(substitute(n)), deparse(substitute(m)), deparse(substitute(u)))
	ng <- length(n)									# 群の数
	stopifnot(ng > 1, length(m) == ng, length(u) == ng, n > 1, floor(n) == n, u > 0)
	if (var.equal) {								# 分散が等しいと仮定する場合
		method <- "二次データに基づく一元配置分散分析（分散が等しいと仮定する場合）"
		nc <- sum(n)								# 全サンプルサイズ
		sw <- sum(u*(n-1))							# 群内変動
		sb <- sum(n*(m-sum(n*m)/nc)^2)						# 群間変動
		ss <- c(sb, sw, sb+sw)							# 平方和（変動）
		df <- c(ng-1, nc-ng, nc-1)						# 自由度
		ms <- ss/df								# 平均平方
		f <- p <- rep(NA, 3)
		f[1] <- ms[1]/ms[2]							# F 値
		p[1] <- pf(f[1], df[1], df[2], lower.tail = FALSE)			# P 値
		anova.table <- cbind(ss, df, ms, f, p)					# 分散分析表
		colnames(anova.table) <- c("平方和", "自由度", "平均平方", "F 値", "P 値")
		rownames(anova.table) <- c("群間", "群内", "合計")
		return(structure(list(statistic=c(F=f[1]), parameter=c("num df"=df[1], "denom df"=df[2]),
			p.value=p[1], method=method, data.name=data.name, anova.table=anova.table),
			class=c("htest", "my.oneway.ANOVA")))
	}
	else {										# 分散が等しいと仮定しない場合（R ではこちらがデフォルト）
		method <- "二次データに基づく一元配置分散分析（分散が等しいと仮定しない場合）"
		w <- n/u
		sum.w <- sum(w)
		m0 <- sum(w*m)/sum.w
		temp <- sum((1-w/sum.w)^2/(n-1))/(ng^2-1)
		f <- sum(w*(m-m0)^2)/((ng-1)*(1+2*(ng-2)*temp))				# F 値
		df1 <- ng-1
		df2 <- 1/(3*temp)
		p <- pf(f, df1, df2, lower.tail=FALSE)				# P 値
		return(structure(list(statistic=c(F=f), parameter=c("num df"=df1, "denom df"=df2), p.value=p,
			method=method, data.name=data.name), class="htest"))
	}
}
# summary メソッド（分散分析表の表示。ただし，等分散を仮定するときのみ）
summary.my.oneway.ANOVA <- function (	obj,
					digits=5)
{
	print(obj$anova.table, na.print="", digits=digits)
}
# ステム・アンド・リーフ
my.stem <- function(	d,
			f=-99)						# f は小数点の移動（元の値を10^f倍した整数部をstemにする）
{
	get.factor <- function(x, minx, maxx)
	{
		for (i in c(1:10, -1:-10)) {
			ll <- trunc(maxx*10^i)-trunc(minx*10^i)
			if (ll >= 2 && ll < 19) {
				 return(10^i)
			}
		}
		return(1)
	}
	DUMMY <- 99							# ダミーの leaf
	MINUS <- -0.1							# -0.xxx などの stem
	d <- d[!is.na(d)]
	f <- ifelse(f == -99, get.factor(d, min(d), max(d)), 10^f)
	temp <- trunc(d*f*10)
	stem <- trunc(temp/10)
	leaf <- abs(temp)-abs(stem*10)
	stem <- ifelse(stem == 0, ifelse(d > 0, 0, MINUS), stem)

	# 跳んでいる stem を補間する
	min.stem <- min(stem)
	max.stem <- max(stem)
	stem <- c(stem, min.stem:max.stem)
	leaf <- c(leaf, rep(DUMMY, max.stem-min.stem+1))

	# -0.xxx と +0.yyy が存在しうるとき
	if (max.stem > 0 && min.stem < 0) {
		stem <- c(stem, MINUS)
		leaf <- c(leaf, DUMMY)
	}

	res <- table(stem, leaf)
	sapply(1:nrow(res),
		function(i) {
			stem.temp <- dimnames(res)$stem[i]
			cat(ifelse(as.numeric(stem.temp) == MINUS, "-0", stem.temp), "| ")
			sapply(1:ncol(res),
				function(le) {
					if (dimnames(res)$leaf[le] != DUMMY) {
						sapply(rep(dimnames(res)$leaf[le], res[i, le]), cat)
					}
				}
			)
			cat("\n")
		}
	)
	cat("stem * ", 1/f, "\n")
}
# 二次データから，二群の平均値の差の検定を行う
my.t.test <- function(	nx,							# 第一群のデータ個数
			mx,							# 第一群の平均値
			ux,							# 第一群の不偏分散
			ny,							# 第二群のデータ個数
			my,							# 第二群の平均値
			uy,							# 第二群の不偏分散
			var.equal = FALSE)					# 等分散を仮定するか
{
	data.name <- sprintf("\nn1 = %s, mean1 = %s, variance1 = %s\nn2 = %s, mean2 = %s, variance2 = %s",
				nx, mx, ux, ny, my, uy)
	if (var.equal) {							# 等分散を仮定するとき，
		method <- "等分散を仮定した，二群の平均値の差の検定"
		df <- nx+ny-2							# 自由度
		v <- ((nx-1)*ux+(ny-1)*uy)/df					# プールした不偏分散
		tstat <- abs(mx-my)/sqrt(v*(1/nx+1/ny))				# 検定統計量
	}
	else {									# 等分散を仮定しないとき
		method <- "ウエルチの方法による，二群の平均値の差の検定"
		tstat <- abs(mx-my)/sqrt(ux/nx+uy/ny)				# 検定統計量
		df <- (ux/nx+uy/ny)^2/((ux/nx)^2/(nx-1)+(uy/ny)^2/(ny-1))	# 自由度（小数点つき）
	}
	P <- 2*pt(tstat, df, lower.tail=FALSE)					# P 値
	names(tstat) <- "t"
	names(df) <- "df"
	return(structure(list(statistic=tstat, parameter=df, p.value=P,
		method=method, data.name=data.name), class="htest"))
}
# 二次データから，二群の等分散性の検定を行う
my.var.test <- function(nx,				# 第一群のデータ個数
			vx,				# 第一群の不偏分散
			ny,				# 第二群のデータ個数
			vy)				# 第二群の不偏分散
{
	data.name <- sprintf("n1 = %s, variance1 = %s, n2 = %s, variance2 = %s", nx, vx, ny, vy)
	method <- "二次データから，二群の等分散性の検定"
	if (vx > vy) {					# f が 1 より大きくなるように
		f <- vx/vy
		df1 <- nx-1
		df2 <- ny-1
	}
	else {
		f <- vy/vx
		df1 <- ny-1
		df2 <- nx-1
	}
	p <- 2*pf(f, df1, df2, lower.tail=FALSE)	# P 値
	dfs <- c("num df"=df1, "denom df"=df2)
	names(f) <- "F"
	return(structure(list(statistic=f, parameter=dfs, p.value=p,
		method=method, data.name=data.name), class="htest"))
}
#####
#
# 相関係数行列を作成し，無相関検定の結果と共に表示する
#
#####


mycor <- function(i,								# 相関係数行列を求める変数が入っているデータフレーム上の列番号または変数名のベクトル
		  df,								# データフレーム
		  use=c("all.obs", "complete.obs", "pairwise.complete.obs"),	# 欠損値の処理法
		  method=c("pearson", "kendall", "spearman"),			# 計算する相関係数の種類
		  latex=TRUE,							# LaTeX 形式で度数分布表を出力する（デフォルトは LaTeX 形式）
		  caption=NULL,							# LaTeX 出力のときの表題
		  label=NULL,							# LaTeX 出力のときのラベル
		  digits=3,							# 相関係数を表示するときの小数点以下の桁数
		  output="",							# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
		  encoding=getOption("encoding"))				# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

# 下請け関数
	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

	outfunc <- function(r, n) {						# 相関係数行列に基づき，若干の計算を行い，結果を出力する
		nr <- nrow(r)							# 相関係数行列のサイズ
		cname <- colnames(r)						# 列名（変数名）をベクトルに入れておく
		rname <- rownames(r)						# 行名（変数名）をベクトルに入れておく（正方行列を扱うので，colnames と同じでよい）
		if (!is.array(n)) {						# n が行列でない（欠損値をリストワイズ除去した場合）
			nall <- sprintf(ifelse(latex, "，$n$=%i", "，n=%i"), n)	# n を nall に文字列としてセットする（latex=TRUE のときは，$n$ とする）
			n <- matrix(n, nr, nr)					# n を 行列にする
		}
		else {
			nall <- ""						# nall を空にセットする
		}
		p <- matrix(0, nr, nr)						# P 値を入れるための行列を作る
		if (method.name == "ケンドールの順位相関係数") {		# ケンドールの順位相関係数の場合の P 値の計算
			z <- abs(r)/sqrt((4*n+10)/(9*n*(n-1)))			# 検定統計量（対称行列だが，無駄に 2 倍の計算をしている）
			p[upper.tri(p, diag=FALSE)] <-				# P 値を計算（こちらは上三角行列だけ）
			    pnorm(z[upper.tri(z, diag=FALSE)], lower.tail=FALSE)*2
		}
		else {								# ピアソンの積率相関係数，スピアマンの順位相関係数の場合の P 値の計算
			t <- abs(r)*sqrt(n-2)/sqrt(1-r^2)			# 検定統計量（対角成分は Inf になるが，P 値を計算するのには何の支障もない）	
			p[upper.tri(p, diag=FALSE)] <-				# P 値を計算（こちらは上三角行列だけ）
			    pt(t[upper.tri(t, diag=FALSE)], n[upper.tri(n, diag=FALSE)]-2, lower.tail=FALSE)*2
		}
		p <- p+t(p)							# 下側三角行列を補完し
		diag(p) <- NA							# 対角成分は NA を代入する
		if (latex) {							# LaTeX 形式で集計結果を出力する
			cat("\n\\begin{table}[htbp]\n", file=output)		# \begin{table}[htbp]
			if (is.null(caption)) {
				cat(sprintf("\\caption{%s行列（欠損値は%s%s）}\n",	# \caption{スピアマンの順位相関係数行列（欠損値はリストワイズ除去，n = xx）}
			    	            method.name, use.name, nall), file=output)	#        --- のようになる
			}
			else {
				cat(sprintf("\\caption{%s}\n", caption), file=output)	# 指定した表題
			}
			if (!is.null(label)) {
				cat(sprintf("\\label{%s}\n", label), file=output)	# 指定したラベル
			}
			cat("\\centering\n", file=output)			# \centering
			cat("\\begin{tabular}{l", rep("c", nr), "} \\hline\n",	# \begin{tabular}{cc…c} \hline
			    sep="", file=output)
			cat("", sprintf("(%03i)", 1:nr), sep=" & ", file=output)# 列名は書けないので (001), (002) … のようにする
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			for (i in 1:nr) {					# 各行ごとに
				cat(sprintf("(%03i)%s",				# (00i)変数名 & 相関係数i1 & 相関係数i2 & …
				    i, rname[i]), sprintf(ifelse(r[i,] < 0, formatl, formatw), r[i,]), sep=" & ", file=output)
				cat("\\\\\n", file=output)			# \\
				cat("$P$値", sprintf(formatw, p[i,]),		# P 値 & 数値i1 & 数値i2 & …
				    sep=" & ", file=output)
				if (use.name == "ペアワイズ除去") {		# 欠損値をペアワイズ除去したのならば，
					cat("\\\\\n", file=output)		# \\
					cat("$n$", n[i,], sep=" & ",file=output)# n & 数値i1 & 数値i2 & …
				}
				cat("\\\\ \\hline\n", file=output)		# \\ \hline
			}
			cat("\\end{tabular}\n", file=output)			# \end{tabular}
			cat("\\end{table}\n", file=output)			# \end{table}
		}
		else {								# tab で区切って出力する
			cat("表　", method.name,				# 表　ピアソンの積率相関係数行列（欠損値はペアワイズ除去）
			    "行列（欠損値は", use.name, nall, "）\n\n",		#        --- のようになる
			    sep="", file=output)
			cat(" ", sprintf("(%03i)", 1:nr),			# 列名は書けないので (001), (002) … のようにする
			    sep="\t", file=output, fill=TRUE)
			for (i in 1:nr) {					# 各行ごとに
				cat(sprintf("(%03i)%s",				# (00i)変数名 & 相関係数i1 & 相関係数i2 & …
				    i, rname[i]), sprintf(formatw, r[i,]),
				    sep="\t", file=output, fill=TRUE)
				cat("P 値", sprintf(formatw, p[i,]),		# P 値 & 数値i1 & 数値i2 & …
				    sep="\t", file=output, fill=TRUE)
				if (use.name == "ペアワイズ除去") {		# 欠損値をペアワイズ除去したのならば，
					cat("n", n[i,], sep="\t",		# n & 数値i1 & 数値i2 & …
					    file=output, fill=TRUE)
				}
			}
		}
	}

# 関数本体
	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

	formatw <- paste("%.", digits, "f", sep="")				# 相関係数の小数点以下の桁数を設定する(Word)
	formatl <- paste("$%.", digits, "f$", sep="")				# 相関係数の小数点以下の桁数を設定する(LaTeX)
	method <- match.arg(method)
	method.name <- switch(match.arg(method),				# 引数 method と変数 method.name の対応付け
			      pearson = "ピアソンの積率相関係数",
			      kendall = "ケンドールの順位相関係数",
			      spearman = "スピアマンの順位相関係数")
	if (method == "pearson" && !all(sapply(i, function(k) is.numeric(df[,k])))) {
		warning("ピアソンの積率相関係数を計算するためには，全ての変数は間隔尺度・比尺度変数であることが必要です")
	}
	else if (!all(sapply(i, function(k) is.numeric(df[,k]) || is.ordered(df[,k])))) {
		warning("順位相関係数を計算するためには，全ての変数は順序尺度・間隔尺度・比尺度変数であることが必要です")
	}
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (match.arg(use) != "pairwise.complete.obs") {			# ペアワイズ除去でない場合（リストワイズ除去の場合）
		use.name <- "リストワイズ除去"					# 欠損値の処理法の名前
		df2 <- df[, i]							# データフレームから分析対象変数を抽出
		df2 <- subset(df2, complete.cases(df2))				# 欠損値 NA を一つも含まないケースを抽出
		n <- nrow(df2)							# サンプルサイズ
		if (n == 0) {
			stop("欠損値をリストワイズ除去した結果，有効なデータがなくなりました。分析に使用する変数を再確認するか，ペアワイズ欠損値処理を選択してやり直してみてください。")
		}
		r <- cor(df2, method=method)					# cor 関数により相関係数行列を求める
		outfunc(r, n)							# 出力関数 outfunc を呼び出す
	}
	else {									# ペアワイズ除去の場合
		use.name <- "ペアワイズ除去"					# 欠損値の処理法の名前
		ln <- length(i)							# 分析対象とする変数の個数（相関係数行列のサイズ）
		r <- n <- matrix(0, ln, ln)
		for (i1 in 1:ln) {
			for (j1 in i1:ln) {					# 全ての変数の組み合わせについて，
				df2 <- df[, c(i[i1], i[j1])]			# データフレームから分析対象となる 2 変数を抽出
				df2 <- subset(df2, complete.cases(df2))		# 2 変数ともに欠損値でないケースを抽出
				n[i1, j1] <- nrow(df2)				# サンプルサイズ
				r[i1, j1] <- cor(df2[1], df2[2], method=method)	# cor 関数により 2 変数間の相関係数を求める
			}
		}
		r <- r+t(r)							# もう半分の三角行列を完成させる
		diag(r) <- 1							# 対角成分は 1 に決まっている
		n <- n+t(n)							# もう半分の三角行列を完成させる
		diag(n) <- diag(n)/2						# このようにして作ると，対角成分は半分にしないといけない
		rownames(r) <- colnames(r) <- colnames(df)[i]			# 相関係数行列の行と列の名前を付ける
		outfunc(r, n)							# 出力関数 outfunc を呼び出す
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
# ニュートン・ラフソン法により一変数間数 f(x)=0 の解を求める
newton <- function(	fun,			# 関数定義
			x1,			# 初期値
			delta=1e-5,		# 数値微分のときの微小数値
			epsilon=1e-14,		# 許容誤差
			max.rotation=100)	# 収束計算上限回数
{
	fun2 <- function(x)			# 数値微分を行う関数
	{
		(fun(x+delta)-fun(x))/delta	# 接線の傾きの近似値
	}
	for (i in 1:max.rotation) {
		x2 <- x1 -fun(x1)/fun2(x1)	# x2 は x1 が改善された解の近似値
		if (abs((x2-x1)/x2) < epsilon) {
			break
		}
		x1 <- x2
	}
	if (i > max.rotation) {
		warning("収束しませんでした")
	}
	x2
}
# 正規分布への適合度の検定
normaldist <- function(	x,						# 度数ベクトル
			b,						# 最初の階級の下限値
			w,						# 階級幅
			a)						# 測定精度
{
	data.name <- deparse(substitute(x))
	method <- "正規分布への適合度の検定"	
	n <- sum(x)							# データ数
	x <- c(0, x, 0)							# 上下にそれぞれ 1 階級を追加する
	k <- length(x)							# 階級数
	mid <- seq(b-w/2, b+k*w-w, w)-a/2				# 級中心

	xbar <- sum(mid*x)/n						# 平均値
	variance <- sum(x*(mid-xbar)^2)/n				# 分散（不偏分散ではない）
	SD <- sqrt(variance)						# 標準偏差
	result <- c("n"=n, "Mean"=xbar, "Variance"=variance, "S.D."=SD)

	z <- ((mid+w/2)-xbar)/SD					# 級限界の標準化得点
	p <- pnorm(z)							# 累積確率
	p[k] <- 1							# 最後の累積確率は 1
	p <- p-c(0, p[-k])						# 各階級の確率
	expectation <- n*p						# 各階級の期待値
	table <- data.frame(mid, x, z, p, expectation)			# 結果をデータフレームにする
	rownames(table) <- paste("c-", 1:k, sep="")

	while (expectation[1] < 1) {					# 期待値が 1 未満の階級を併合
		x[2] <- x[2]+x[1]
		expectation[2] <- expectation[2]+expectation[1]
		x <- x[-1]
		expectation <- expectation[-1]
		k <- k-1
	}
	while (expectation[k] < 1) {					# 期待値が 1 未満の階級を併合
		x[k-1] <- x[k-1]+x[k]
		expectation[k-1] <- expectation[k-1]+expectation[k]
		x <- x[-k]
		expectation <- expectation[-k]
		k <- k-1
	}
	chisq <- sum((x-expectation)^2/expectation)			# カイ二乗統計量
	df <- k-3							# 自由度
	p <- pchisq(chisq, df, lower.tail=FALSE)			# P 値
	names(chisq) <- "X-squared"
	names(df) <- "df"
	return(structure(list(statistic=chisq, parameter=df, p.value=p,
		estimate=c(n=n, Mean=xbar, Variance=variance, S.D.=SD),
		method=method, data.name=data.name, table=table),
		class=c("htest", "normaldist")))
}
# summary メソッド
summary.normaldist <- function(	obj,					# normaldist が返すオブジェクト
				digits=5)				# 表示桁数
{
	table <- obj$table
	colnames(table) <- c("級中心", "度数", "標準化得点", "確率", "期待値")
	cat("\n適合度\n\n")
	print(table, digits=digits, row.names=FALSE)
}
# plot メソッド
plot.normaldist <- function(	obj,					# normaldist が返すオブジェクト
				xlab="", ylab="Frequency",		# 軸の名称
				...)					# plot に渡す引数
{
	d <- obj$table
	class <- d$mid
	w <- diff(class)[1]
	class <- class-w/2
	k <- length(class)
	yo <- d$x
	ye <- d$expectation
	plot(c(class[1], class[k]+w), c(0, max(c(yo, ye))), type="n", xlab=xlab, ylab=ylab, ...)
	rect(class, 0, class+w, yo, col="grey")
	Mean <- obj$estimate[2]
	SD <- obj$estimate[4]
	x <- seq(class[1], class[k]+w, length=2000)
	y <- dnorm(x, mean=Mean, sd=SD)*sum(yo)*w
	lines(x, y)
	x <- d$mid
	y <- dnorm(x, mean=Mean, sd=SD)*sum(yo)*w
	points(x, y, pch=3)
}
# 正規確率紙に累積相対度数をプロットする
npp <- function(y,					# 度数ベクトル
		x=NULL,					# 階級代表値
		plt=TRUE,				# データ点をプロットする
		xlab=NULL,				# 横軸ラベル
		ylab=NULL,				# 縦軸ラベル
		main=NULL)				# 図のタイトル
{
	if (length(y) < 3) {				# 階級数が2以下のときには正規確率紙のみを出力する
		y <- 1:11
		plt <- FALSE
	}
	y <- cumsum(y)/sum(y)				# 累積相対度数
	if (is.null(x)) {
		x = seq(along=y)
	}
	if (is.null(xlab)) xlab <- "観察値"
	if (is.null(ylab)) ylab <- "累積相対度数"
	if (is.null(main)) main <- "正規確率紙"
	probs <- c(0.01, 0.1, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 99.9, 99.99)/100
	plot(c(x[1], x[length(x)-1]), qnorm(c(probs[1], probs[17])), type="n", xaxt="n", yaxt="n",
		xlab=xlab, ylab=ylab, main=main)
	abline(h=qnorm(probs), v=x, col="grey")		# 格子線を引く
	if (plt) {					# データ点をプロットする
		points(x, qnorm(y), type="b")
		text(x, qnorm(y), round(y, digit=3)*100, pos=1)
	}
	axis(1, at=x)					# 横軸を描く
	axis(2, at=qnorm(probs), labels=probs*100)	# 縦軸を描く
}
# 正規確率紙に累積相対度数をプロットする
npp2 <- function(x)		# データベクトル
{
	x <- x[!is.na(x)]	# 欠損値を持つケースを除く
	n <- length(x)		# データの個数
	x <- sort(x)		# 昇順に並べ替える
	y <- (1:n-0.5)/n	# 累積相対度数
	probs <- c(0.01, 0.1, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 99.9, 99.99)/100
	plot(c(x[1], x[n]), qnorm(c(probs[1], probs[17])), type="n", yaxt="n",
		xlab="Observed Value", ylab="Cumulative Percent",
		main="Normal Probability Paper")
	points(x, qnorm(y))
	axis(2, qnorm(probs), probs*100)
}
# オッズ比を計算する
odds.ratio <- function(	a, b, c, d,		# 4 つのセルの観察値
			correct=FALSE)		# 補正をするとき TRUE にする
{
	cl <- function(x)
	{
		or*exp(c(1, -1)*qnorm(x)*sqrt(1/a+1/b+1/c+1/d))
	}
	if (correct || a*b*c*d == 0) {		# a,b,c,d のいずれかが 0 のときには必ず補正する
		a <- a+0.5
		b <- b+0.5
		c <- c+0.5
		d <- d+0.5
	}
	or <- a*d/(b*c)
	conf <- rbind(cl90=cl(0.05), cl95=cl(0.025), cl99=cl(0.005), cl999=cl(0.0005))
	conf <- data.frame(conf)
	colnames(conf) <- paste(c("下側","上側"), "信頼限界値", sep="")
	rownames(conf) <- paste(c(90, 95, 99, 99.9), "%信頼区間", sep="")
	list(or=or, conf=conf)
}
oneway.ANOVA <- function(x, g, verbose=TRUE) { # x: データベクトル，g: 群変数ベクトル
# 基礎統計量
  d <- data.frame(x, g)
  level <- levels(d$g)
  k <- length(level)
  n <- aggregate(x ~ g, d, length)[,2]
  N <- sum(n)
  Mean <- aggregate(x ~ g, d, mean)[,2]
  Var <- aggregate(x ~ g, d, var)[,2]
  SD <- aggregate(x ~ g, d, sd)[,2]
  SE <- SD/sqrt(n)
  result <- data.frame(水準=level, n, Mean,
                       SD, "Mean-SD"=Mean-SD, "Mean+SD"=Mean+SD,
                       SE, "Mean-SE"=Mean-SE, "Mean+SE"=Mean+SE, check.names=FALSE)
# デザインプロット
  if (verbose) {
    x.pos <- seq_along(level)
    plot(x.pos, Mean, xlim=c(0.5, k+0.5), ylim=range(Mean-SD, Mean+SD), pch=16,
         xlab="", xaxt="n", ylab="")
    axis(1, at=x.pos, labels=level)
    arrows(x.pos, Mean-SD, x.pos, Mean+SD, length=0.10, angle=90, code=3, lty=3)
    arrows(x.pos, Mean-SE, x.pos, Mean+SE, length=0.05, angle=90, code=3, lwd=2)
  }
# 等分散性の検定
  r.bartlett <- bartlett.test(d$x, d$g)
  r.Levene <- Levene.test(d$x, d$g)
  r.eqvar <- data.frame(手法=c("Bartlett", "Levene"), 分布=c("chi square", "F"),
               統計量=c(r.bartlett$statistic, r.Levene$statistic),
               自由度=c(r.bartlett$parameter, 
                        sprintf("(%i, %i)", r.Levene$parameter[1], r.Levene$parameter[2])),
               "P 値"=c(r.bartlett$p.value, r.Levene$p.value))
# 分散分析（等分散性を仮定する）
  r.oneway <- oneway.test(x ~ g, d, var.equal=TRUE)
  gm <- mean(d$x)
  St <- sum((d$x-gm)^2)
  Sb <- sum(n*(Mean-gm)^2)
  Sw <- St-Sb
  SS <- c(Sb, Sw, St)
  df <- c(k-1, N-k, N-1)
  MS <- SS/df
  F <- c(MS[1]/MS[2], NA, NA)
  p <- c(pf(F[1], df[1], df[2], lower.tail=FALSE), NA, NA)
  MS[3] <- NA
  r.ANOVA <- data.frame(要因=c("級間", "級内", "全体"), 平方和=c(Sb, Sw, St),
                        自由度=df, 平均平方=MS,  "F 値"=F, "P 値"=p, check.names=FALSE)
# 分散分析（等分散性を仮定しない）
  r.Welch <- oneway.test(x ~ g, d)
  r.BF <- Brown.Forsythe.test(x, g)
  r.ANOVA2 <- data.frame(手法=c("Welch", "Brown-Forsythe"),
                         "F 値"=c(r.Welch$statistic, r.BF$statistic),
                         "第 1 自由度"=c(k-1, r.BF$parameter[1]),
                         "第 2 自由度"=c(r.Welch$parameter[2], r.BF$parameter[2]),
                         "P 値"=c(r.Welch$p.value, r.BF$p.value), check.names=FALSE)
# 多重比較
  r.tukey <- tukey(d$x, d$g)$Tukey
  lc <- t(sapply(rownames(r.tukey), function(s) as.integer(unlist(strsplit(s, ":")))))
  diff <- Mean[lc[,2]]-Mean[lc[,1]]
  r.GH <- tukey(d$x, d$g, method="G")$Games.Howell
  r.multi <- data.frame(群1=level[lc[,1]], 群2=level[lc[,2]],
                        平均1=Mean[lc[,1]], 平均2=Mean[lc[,2]], 差=diff,
                        差の標準誤差=diff/r.tukey[,1], "t 値(Tukey)"=r.tukey[,1],
                        "P 値"=r.tukey[,2], "t 値(Games-Howell)"=r.GH[,1],
                        "自由度"=r.GH[,2], "P 値"=r.GH[,3], check.names=FALSE)
  if (verbose) {
    cat("\n基礎統計量\n\n")
    print(result, row.names=FALSE)
    cat("\n等分散性の検定\n\n")
    print(r.eqvar, row.names=FALSE)
    cat("\n分散分析（等分散性を仮定する）\n\n")
    print(r.ANOVA, row.names=FALSE)
    cat("\n分散分析（等分散性を仮定しない）\n\n")
    print(r.ANOVA2, row.names=FALSE)
    cat("\n多重比較\n\n")
    print(r.multi, row.names=FALSE)
  }
  invisible(list(result, r.eqvar, r.ANOVA, r.ANOVA2, r.multi))
}
# 二本の直線による，折れ線回帰を行う
oresen <- function(	x,			# 独立変数
			y)			# 従属変数
{
	ss <- function(par)			# 暫定条件下の残差平方和を求める
	{
		a <- par[1]			# 交点の x 座標
		b <- par[2]			# 交点の y 座標
		c <- par[3]			# 左側の直線の傾き
		d <- par[4]			# 右側の直線の傾き
		xl <- x[x < a]			# データの左部分
		xr <- x[x >= a]			# データの右部分
		if (length(xl) == 0 ||          # 左右いずれかにデータがないときには無限大を返す
		    length(xr) == 0) {
			return(Inf)
		}
		yl <- y[x < a]
		yr <- y[x >= a]
		yle <- c*(xl-a)+b		# 左部分の予測値
		yre <- d*(xr-a)+b		# 左部分の予測値
		retv <- sum((yl-yle)^2)+
			sum((yr-yre)^2)		# 残差平方和の和
		return(retv)
	}
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	OK <- complete.cases(x, y)		# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	par <- c(mean(x), mean(y), 1, 1)	# 初期値
	ans <- optim(par, ss, control=list(maxit=1000))
	obj <- list(names.xy=names.xy, x=x, y=y, par=ans$par, residuals=ans$value)
	class(obj) <- "oresen"
	return(obj)
}
# print メソッド
print.oresen <- function(obj)
{
	cat(sprintf("残差平方和 = %g\n", obj$residuals))
	par <- obj$par
	cat(sprintf("交点座標 = ( %g, %g )\n", par[1], par[2]))
	cat(sprintf("切片 = %g, 傾き = %g\n", -par[3]*par[1]+par[2], par[3])) # 左側の回帰直線の式
	cat(sprintf("切片 = %g, 傾き = %g\n", -par[4]*par[1]+par[2], par[4])) # 右側の回帰直線の式
}
# plot メソッド
plot.oresen <- function(obj,                    # oresen オブジェクト
			xlab=obj$names.xy[1],	# x 軸の名前
			ylab=obj$names.xy[2],	# y 軸の名前
			col1="red",		# 左側のデータ点と回帰直線の描画色
			col2="blue",		# 右側のデータ点と回帰直線の描画色
			...)			# plot 関数などに渡す引数
{
	par <- obj$par
	x <- obj$x
	y <- obj$y
	plot(x,  y, xlab=xlab, ylab=ylab, col=ifelse(x < par[1], col1, col2), ...)
	abline(-par[3]*par[1]+par[2], par[3], col=col1, ...)
	abline(-par[4]*par[1]+par[2], par[4], col=col2, ...)
}
# ライアンの方法とチューキーの方法による比率の対比較
p.multi.comp <- function(	n,				# 標本サイズ
				r,				# 陽性サイズ
				alpha=0.05,			# 有意水準
				method=c("ryan", "tukey"))	# 方法
{
	printf <- function(fmt, ...)
	{
		cat(sprintf(fmt, ...))
	}

	check <- function(s, b)					# 検定しようとしている二群が，それまでに有意でないとされた二群に挟まれているか
	{
		if (ns.n > 1) {
			for (i in 1:ns.n) {
				if (ns.s[i] <= s && s <= ns.b[i] && ns.s[i] <= b && b <= ns.b[i]) {
					return(FALSE)		# 検定するまでもなく有意でないとする
				}
			}
		}
		return(TRUE)					# 検定しなくてはならない
	}

	k <- length(n)						# 群の数
	stopifnot(k == length(r), k == length(n), n > 0, r >= 0, r <= n, floor(n) == n, floor(r) == r)
	method <- match.arg(method)				# 引数の補完
	o <- order(r/n)						# 標本比率の大きさの順位
	sr <- r[o]
	sn <- n[o]
	sr.sn <- sr/sn						# 割合
	num.significant <- ns.n <- 0
	ns.s <- ns.b <- numeric(k*(k-1)/2)			# 有意でない群の記録用
	for (m in k:2) {
		for (small in 1:(k-m+1)) {
			big <- small+m-1
			if (check(small, big)) {
				prop <- sum(sr[small:big]) / sum(sn[small:big])		# 推定母比率
				se <- sqrt(prop*(1-prop)*(1/sn[small]+1/sn[big]))	# 標準誤差
				diff <- sr.sn[big]-sr.sn[small]				# 比率の差
				if (method == "ryan") { # Ryan の方法
					nominal.alpha <- 2*alpha/(k*(m-1))		# 名義的有意水準
					rd <- se*qnorm(nominal.alpha/2, lower.tail=FALSE)	# 差があると見なせる差の大きさ
					z <- if (se == 0) 0 else diff/se
					p <- pnorm(z, lower.tail=FALSE)*2
					result <- p <= nominal.alpha
				}
				else { # Tukey の方法
					if (m == k) {		# 最大群と最小群の比較のとき
						qq <- q <- qtukey(alpha, k, Inf, lower.tail=FALSE)
					}
					else {                  # 上記以外の群の比較のとき
						qq <- (q+qtukey(alpha, m, Inf, lower.tail=FALSE))/2
					}
					WSD <- se*qq/sqrt(2)
					result <- diff != 0 && diff >= WSD
				}
				if (result) { # 有意であるとき
					num.significant <- 1
					printf("p[%2i]=%7.5f vs. p[%2i]=%7.5f : diff.= %7.5f, ",
						o[small], sr.sn[small], o[big], sr.sn[big], diff)
					if (method == "ryan") {
						printf("RD=%7.5f : P=%7.5f, alpha'=%7.5f\n", rd, p, nominal.alpha)
					}
					else {
						printf("WSD=%7.5f\n", WSD)
					}
				}
				else {				# 有意でないとき
					ns.n <- ns.n+1
					ns.s[ns.n] <- small
					ns.b[ns.n] <- big
				}
			}
		}
	}
	if (num.significant == 0) {				# 有意差のある群は一つもなかった
		print("Not significant at all.")
	}
}
# 対応のある平均値の差の検定
paired.t.test <- function(	x, y,			# 一次データから検定を行うときは対応のあるベクトル，
          						# 二次データから検定を行うときは対応のある二つの平均値
				ux2=NULL, uy2=NULL,	# 二次データから検定を行うときは対応のある二つの標準偏差
				r=NULL,			# 二次データから検定を行うとき，対応のあるデータ間の相関係数
				n=NULL)			# 二次データから検定を行うとき，データの組数
{
	method <- "対応のある平均値の差の検定"
	if (is.null(ux2)) {				# 一次データについて検定
		data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
		OK <- complete.cases(x, y)		# 欠損値を持つケースを除く
		x <- x[OK]
		y <- y[OK]
		n <- length(x)				# サンプルサイズ
		x <- x-y				# 差を取って，
		t <- abs(mean(x))/(sd(x)/sqrt(n))	# 検定
	}
	else {						# 二次データについて検定
		data.name <- sprintf("\nmean1 = %s, variance1 = %s\nmean2 = %s, variance2 = %s\nr = %s, n = %s",
					x, ux2, y, uy2, r, n)
		t <- abs(x-y)/sqrt((ux2+uy2-2*r*sqrt(ux2*uy2))/n)
	}
	df <- n-1					# 自由度
	p <- pt(t, df, lower.tail=FALSE)*2		# P 値
	return(structure(list(statistic=c(t=t), parameter=c(df=df), p.value=p,
		method=method, data.name=data.name), class="htest"))
}
# 比率の差（分布の差）の多重比較（prop.test, fisher.test, chisq.test を用いる）
pairwise.prop2.test <- function (	x,					# 比率の差（2 カテゴリー）ならば，ベクトルまたは列数 2 の行列。または，列数 2 以上の行列
					n,					# 比率の差（2 カテゴリー）ならば，ベクトル。それ以外の場合は無視される
					p.adjust.method = p.adjust.methods,	# P 値の調整法
					test.function=prop.test, ...) 		# 下請けに使う関数名 prop.test の他に，fisher.test, chisq.test が使える
{
    p.adjust.method <- match.arg(p.adjust.method)
    METHOD <- deparse(substitute(test.function))
    DNAME <- deparse(substitute(x))
    if (is.matrix(x)) {
        if (ncol(x) < 2) 
            stop("'x' must have at least 2 columns")
    } else if (is.vector(x) && is.vector(n))
    	x <- cbind(x, n-x)
    else
    	stop("'x' must be a matrix, or 'x', and 'n' must be a vector")
    if (nrow(x) < 2) 
        stop("too few groups")
    compare.levels <- function(i, j) {
        test.function(x[c(i, j),], ...)$p.value					# test.function で使用する検定を使い分ける
    }
    level.names <- names(x)
    if (is.null(level.names)) 
        level.names <- seq_along(1:nrow(x))
    PVAL <- pairwise.table(compare.levels, level.names, p.adjust.method)	# R の pairwise.*.test は有意水準ではなく P 値を調整するので解釈時に注意
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
        p.adjust.method = p.adjust.method)
    class(ans) <- "pairwise.htest"
    ans
}
# 代表値の差の多重比較（対比較）
pairwise.wilcox.test <- function (x, g, p.adjust.method = p.adjust.methods,
	paired = FALSE, exact=TRUE, ...)
{
	DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
	g <- factor(g)
	p.adjust.method <- match.arg(p.adjust.method)
	METHOD <- if (paired) "paired wicoxon tests" else "wilcoxon tests"
	compare.levels <- function(i, j) {
		xi <- x[as.integer(g) == i]
		xj <- x[as.integer(g) == j]
		if (paired) {
			ok <- complete.cases(xi, xj)
			xi <- xi[ok]
			xj <- xj[ok]
		}
		else {
			xi <- xi[!is.na(xi)]
			xj <- xj[!is.na(xj)]
		}
		if (exact && !paired) {
			library(coin)
			x <- c(xi, xj)
			g <- factor(rep(1:2, c(length(xi), length(xj))))
			pvalue(wilcox_test(x ~ g, distribution="exact"))
		}
		else
			wilcox.test(xi, xj, paired = paired, ...)$p.value
	}
	PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
	ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
		p.adjust.method = p.adjust.method)
	class(ans) <- "pairwise.htest"
	ans
}
# 偏相関係数行列
partial.cor <- function(x)					# データ行列
{
	x <- subset(x, complete.cases(x))			# 欠損値を持つケースを除く
	i <- solve(cor(x))					# 相関係数行列の逆行列
	d <- diag(i)						# 対角成分
	i <- -i/sqrt(outer(d, d))				# 偏相関係数行列
	diag(i) <- NA						# 対角成分は未定義
	rownames(i) <- colnames(i) <- paste("Var", 1:ncol(x))
	return(i)
}
# 一対比較データを双対尺度法で分析する
pc.dual <- function(	F,					# 一対比較データ
			one.two=TRUE,				# 1/2 で入力されているとき 1/-1 に変換する
			col.names=NULL)				# 評価対象名
{
	F <- data.matrix(F)					# データフレームも行列にする
	if (one.two) {						# 1/2 で入力されているとき 1/-1 に変換する
		F[F == 2] <- -1
	}
	N <- nrow(F)						# 被調査者数
	if (is.null(rownames(F))) {				# 被験者の名前が与えられていないとき
		row.names <- paste("Row", 1:N, sep="-")		# 行ラベルの補完
	}
	n <- (1+sqrt(1+8*ncol(F)))/2				# 比較対象とされるものの数
	if (is.null(col.names)) {				# 比較対象の名前が与えられていないとき
		col.names <- paste("Col", 1:n, sep="-")		# 列ラベルの補完
	}
	x <- combn(n, 2)					# モデル行列作成の添え字
	nc <- ncol(x)
	A <- matrix(0, nc, n)					# モデル行列作成
	A[cbind(1:nc, x[1,])] <- 1
	A[cbind(1:nc, x[2,])] <- -1
	E <- F%*%A
	Hn <- t(E)%*%E/(N*n*(n-1)^2)
	ans <- eigen(Hn)					# 固有値・固有ベクトルを求める
	ne <- nrow(Hn)-1					# 有効な固有値・固有ベクトルの個数
	eta2 <- ans$values[1:ne]				# 固有値（相関比の二乗）
	eta <- sqrt(eta2)					# 相関比
	contribution <- eta2/sum(ans$values[1:ne])*100		# 寄与率
	cumcont <- cumsum(contribution)				# 累積寄与率
	result <- rbind(eta2, eta, contribution, cumcont)	# 結果
	dimnames(result) <- list(c("eta square", "correlation", "contribution", "cumulative contribution"),
				   paste("Axis", 1:ne, sep="-"))
	W <- ans$vectors[, 1:ne, drop=FALSE]			# 固有ベクトル
	col.score <- W*sqrt(n)					# 列スコア
	col.score2 <- t(t(col.score)*eta)			# 相関比で重み付けした列スコア
	row.score2 <- t(t(E%*%W/sqrt(n)/(n-1)))			# 相関比で重み付けした行スコア
	row.score <- t(t(row.score2)/eta)			# 行スコア
	colnames(col.score) <- colnames(row.score) <- colnames(result)
	rownames(col.score) <- col.names
	rownames(row.score) <- row.names
	dimnames(col.score2) <- dimnames(col.score)
	dimnames(row.score2) <- dimnames(row.score)
	result <- list(	result=result,
			row.score=row.score,
			col.score=col.score, 
			row.score.weighted=row.score2, 
			col.score.weighted=col.score2)
	class(result) <- "dual"					# summary, plot メソッドがある
	invisible(result)
}
# 主成分分析
pca <- function(dat)					# データ行列
{
	if (is.null(rownames(dat))) rownames(dat) <- paste("#", 1:nrow(dat), sep="")
	dat <- subset(dat, complete.cases(dat))		# 欠損値を持つケースを除く
	nr <- nrow(dat)					# サンプルサイズ
	nc <- ncol(dat)					# 変数の個数
	if (is.null(colnames(dat))) {
		colnames(dat) <- paste("X", 1:nc, sep="")
	}
	vname <- colnames(dat)
	heikin <- colMeans(dat)				# 各変数の平均値
	bunsan <- apply(dat, 2, var)			# 各変数の不偏分散
	sd <- sqrt(bunsan)				# 各変数の標準偏差
	r <-cor(dat)					# 相関係数行列
	result <- eigen(r)				# 固有値・固有ベクトルを求める
	eval <- result$values				# 固有値
	evec <- result$vectors				# 固有ベクトル
	contr <- eval/nc*100				# 寄与率（％）
	cum.contr <- cumsum(contr)			# 累積寄与率（％）
	fl <- t(sqrt(eval)*t(evec))			# 主成分負荷量
	fs <- scale(dat)%*%evec*sqrt(nr/(nr-1))		# 主成分得点
	names(heikin) <- names(bunsan) <- names(sd) <-
		rownames(r) <- colnames(r) <- rownames(fl) <- colnames(dat)
	names(eval) <- names(contr) <- names(cum.contr) <-
		colnames(fl) <- colnames(fs) <- paste("PC", 1:nc, sep="")
	return(structure(list(mean=heikin, variance=bunsan,
			standard.deviation=sd, r=r,
			factor.loadings=fl, eval=eval,
			contribution=contr,
			cum.contribution=cum.contr, fs=fs), class="pca"))
}
# print メソッド
print.pca <- function(	obj,				# pca が返すオブジェクト
			npca=NULL,			# 表示する主成分数
			digits=3)			# 結果の表示桁数
{
	eval <- obj$eval
	nv <- length(eval)
	if (is.null(npca)) {
		npca <- sum(eval >= 1)
	}
	eval <- eval[1:npca]
	cont <- eval/nv
	cumc <- cumsum(cont)
	fl <- obj$factor.loadings[, 1:npca, drop=FALSE]
	rcum <- rowSums(fl^2)
	vname <- rownames(fl)
	max.char <- max(nchar(vname), 12)
	fmt1 <- sprintf("%%%is", max.char)
	fmt2 <- sprintf("%%%is", digits+5)
	fmt3 <- sprintf("%%%i.%if", digits+5, digits)
	cat("\n主成分分析の結果\n\n")
	cat(sprintf(fmt1, ""),
	    sprintf(fmt2, c(sprintf("PC%i", 1:npca), "  Contribution")), "\n", sep="", collapse="")
	for (i in 1:nv) {
		cat(sprintf(fmt1, vname[i]),
		    sprintf(fmt3, c(fl[i, 1:npca], rcum[i])),
		    "\n", sep="", collapse="")
	}
	cat(sprintf(fmt1, "Eigenvalue"),   sprintf(fmt3, eval[1:npca]), "\n", sep="", collapse="")
	cat(sprintf(fmt1, "Contribution"), sprintf(fmt3, cont[1:npca]), "\n", sep="", collapse="")
	cat(sprintf(fmt1, "Cum.contrib."), sprintf(fmt3, cumc[1:npca]), "\n", sep="", collapse="")
	
}
# summary メソッド
summary.pca <- function(obj,				# pca が返すオブジェクト
			digits=5)			# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.pca <- function(	obj,				# pca が返すオブジェクト
			which=c("loadings", "scores"),	# 主成分負荷量か主成分得点か
			pc.no=c(1,2),			# 描画する主成分番号
			ax=TRUE,			# 座標軸を描き込むかどうか
			label.cex=0.6,			# 主成分負荷量のプロットのラベルのフォントサイズ
			...)				# plot に引き渡す引数
{
	which <- match.arg(which)
	if (which == "loadings") {
		d <- obj$factor.loadings
	}
	else {
		d <- obj$fs
	}
	label <- sprintf("第%i主成分", pc.no)
	plot(d[, pc.no[1]], d[, pc.no[2]], xlab=label[1], ylab=label[2], ...)
	if (which == "loadings") {
		text(d[, pc.no[1]], d[, pc.no[2]], rownames(obj$factor.loadings), pos=1, cex=label.cex)
	}
	abline(h=0, v=0)
}
# 画像を pdf ファイルへ出力するためのデフォルト関数
pdf2 <- function(	fn,					# 出力ファイル名
			of=FALSE)				# 一つのグラフを一つのファイルに出力
{
	pdf(file=fn, width=800/72, height=600/72, onefile=of)	# Mac において，横800ドット×縦600ドットの画像
}
# usage: permutation.oneway.test(list(c(36.7, 52.4, 65.8), c(45.7, 61.9, 65.3), c(52.6, 76.6, 81.3)))
permutation.oneway.test <- function(x)
{
	n <- sapply(x, length)
	g <- factor(rep(1:length(x), n))
	z <- unlist(x)
	library(e1071)
	perm <- permutations(sum(n))
	p0 <- oneway.test(z~g)$p.value
	ps <- apply(perm, 1, function(i) oneway.test(z[i]~g)$p.value)
	invisible(list(p0=p0, permutation.p=mean(ps < p0+1e-10), ps=ps))
}
# 独立2標本の並べ替え検定
permutation.test <- function(	x,		# 第 1 群のデータ
				y,		# 第 2 群のデータ
				FUNC, ...)	# 検定関数および，その引数
{
	x <- x[!is.na(x)]			# 欠損値を除く
	y <- y[!is.na(y)]			# 欠損値を除く
	s <- c( FUNC(x, y, ...)$statistic,	# 観察値の場合の統計量1
		FUNC(y, x, ...)$statistic)	# 観察値の場合の統計量2
	s1 <- min(s)
	s2 <- max(s)
	z <- c(x, y)				# データをプールする
	perm <- combn(length(z), length(x),	# 全ての組み合わせについて
		function(g) {			# 検定を行い，
			r <- FUNC(z[g], z[-g], ...)$statistic # 検定統計量が
			r <= s1 || r >= s2}	# 観察値の場合より極端なら TRUE を返す
		)
	return(mean(perm))			# P 値（= 全ての組み合わせ中での TRUE の数の割合）
}
# 因子分析
pfa <- function(dat,							# データ行列または相関係数行列
		method=c("Varimax", "biQuartimax", "Quartimax", "Equimax", "None"),	# 回転法
		eps1=1e-5,						# 共通性の収束限界
		eps2=1e-5,						# バリマックス回転の収束限界
		max1=999,						# 共通性の推定のための収束計算の上限回数
		max2=999,						# バリマックス回転を行う上限回数
		factors=0)
{
	method <- match.arg(method)					# 回転法を補完
	dat <- subset(dat, complete.cases(dat))				# 欠損値を持つケースを除く
	nr <- nrow(dat)							# ケース数
	nc <- ncol(dat)							# 変数の個数
	if (is.null(colnames(dat))) {					# 変数名が無いときには仮の名前を付ける
		colnames(dat) <- paste("Var", 1:nc, sep=".")
	}
	vnames <- colnames(dat)						# 変数名を記録
	if (nr != nc && is.null(rownames(dat))) {			# ケース名がないときには仮の名前を付ける
		rownames(dat) <- paste("Obs", 1:nr, sep=".")
	}
	cnames <- rownames(dat)						# ケース名を記録
	r0 <- r <- if (nr == nc) dat else cor(dat)			# 与えられたのがデータ行列なら，相関係数行列を計算する
	communality0 <- 1-1/diag(solve(r))				# 共通性の初期値（SMC）
	diag(r) <- communality0						# 相関係数行列の対角成分を共通性で置き換える
	result <- eigen(r)						# 固有値・固有ベクトルを求める
	eval <- result$values						# 固有値
	evec <- result$vectors						# 固有ベクトル
	if (factors == 0) {						# 因子数の指定がないときには，
		factors <- sum(eval >= 1)				# 1 以上の固有値の数とする
	}
	converged <- FALSE						# 共通性の収束計算を行う
	for (i in 1:max1) {						# 上限回数まで繰り返す
		eval <- eval[1:factors]					# 必要な因子数まで，固有値を取る
		evec <- evec[,1:factors]				# 必要な因子数まで
		evec <- t(sqrt(eval)*t(evec))				# 因子負荷量を計算
		r <- r0							# 相関係数行列を復帰
		communality <- rowSums(evec^2)				# 共通性を計算
		if (all(abs(communality-communality0) < eps1)) {	# 共通性の変化が限界以内になれば，
			converged <- TRUE				# 収束したとして，
			break						# ループを抜ける
		}
		communality0 <- communality				# 現在の共通性を保存
		diag(r) <- communality					# 相関係数行列の対角成分を置き換える
		result <- eigen(r)					# 固有値・固有ベクトルを求める
		eval <- result$values					# 固有値
		evec <- result$vectors					# 固有ベクトル
	}
	if (converged == FALSE) {					# 収束フラッグが FALSE なら，
		warning("Not converged.")				# 収束しなかったと警告する
	}
	else if (any(communality >= 1)) {				# 1 を超える共通性があったら，
		warning("Communality >= 1.")				# 警告する
	}

	if (factors == 1 || method == "None") {				# 因子数が 1 であるか，因子軸の回転をしないのなら，
		w <- solve(r0)%*%evec					# 因子得点係数を計算する
		scores <- (scale(dat)*sqrt(nr/(nr-1)))%*%w		# 因子得点を計算する
		rownames(evec) <- names(communality) <- vnames		# 名前を付けて，結果を返す
		rownames(scores) <- cnames
		colnames(scores) <- colnames(evec) <- names(eval) <- paste("FAC", 1:factors, sep=".")
		return(structure(list(rotation=method, correlation.matrix=r0, communality=communality,
			before.rotation.fl=evec, before.rotation.eval=eval, scores=scores), class="pfa"))
	}
	else {								# 因子軸を回転するなら，
		fl <- evec/sqrt(communality)				# 因子負荷量を共通性で基準化
		eig <- numeric(factors)
		ov <- 0
		wt <- switch (method,					# 回転方法によって，重みを変える
			"Varimax" = 1,
			"biQuartimax" = 0.5,
			"Quartimax" = 0,
			"Equimax"	= nc/2)
		fnp <- nc
		for (loop in 1:max2) {					# 回転の上限回数まで収束計算
			for (k in 1:(factors-1)) {
				for (k1 in (k+1):factors) {
					x <- fl[,k]
					y <- fl[,k1]
					xy <- x^2-y^2
					a <- sum(xy)
					b <- 2*sum(x*y)
					c <- sum(xy^2-4*x^2*y^2)
					d <- 4*sum(x*y*xy)
	
					dd <- d-2*wt*a*b/fnp
					theta <- atan(dd/(c-wt*(a^2-b^2)/fnp))
					if(sin(theta)*dd < 0) {
						if (theta > 0) {
							theta <- theta-pi
						}
						else {
							theta <- theta+pi
						}
					}
					theta <- theta/4
					cs <- cos(theta)
					si <- sin(theta)
					fljk <- fl[,k]
					fljk1 <- fl[,k1]
					fl[,k] <- fljk*cs+fljk1*si
					fl[,k1] <- fljk1*cs-fljk*si
				}
			}
			v <- sum((t(fl)^2-colSums(fl^2)*wt/fnp)^2)
	
			if (abs(v-ov) < eps2) {				# 収束したら，
				break					# ループから抜ける
			}
			ov <- v
		}
		fl <- fl*sqrt(communality)				# 因子負荷量
		w <- solve(r0)%*%fl					# 因子得点係数
		scores <- (scale(dat)*sqrt(nr/(nr-1)))%*%w		# 因子得点
		eval2 <- colSums(fl^2)					# 因子負荷量の二乗和
		rownames(evec) <- rownames(fl) <- names(communality) <- vnames
		rownames(scores) <- cnames
		colnames(scores) <- colnames(evec) <- colnames(fl) <- names(eval) <- names(eval2) <- paste("FAC", 1:factors, sep=".")
		return(structure(list(rotation=method, correlation.matrix=r0, communality=communality,
			before.rotation.fl=evec, before.rotation.eval=eval,
			after.rotation.fl=fl, after.rotation.eval=eval2, scores=scores), class="pfa"))
	}
}
# print メソッド
print.pfa <- function(	result,						# pfa が返すオブジェクト
			before=FALSE)
{
	communality <- result$communality
	vnames <- sapply(names(communality), function(i) substring(i, 1, min(nchar(i), 7)))
	if (before || is.null(result$after.rotation.fl)) {
		fl <- result$before.rotation.fl
		eval <- result$before.rotation.eval
		label <- "E-value"
		if (result$rotation == "None") {
			printf("\nResult without rotation\n\n")
		}
		else {
			printf("\nBefore %s rotation\n\n", result$rotation)
		}
	}
	else {
		fl <- result$after.rotation.fl
		eval <- result$after.rotation.eval
		label <- "Sq.Sum"
		printf("\nAfter %s rotation\n\n", result$rotation)
	}
	nv <- nrow(fl)
	nc <- ncol(fl)
	cat("       ", sprintf(" Fac.%02i", 1:nc), " Communality\n", sep="")
	for (i in 1:nv) {
		cat(sprintf("%-7s", vnames[i]), sprintf("%7.3f", fl[i,]), sprintf("%7.3f\n", communality[i]), sep="")
	}
	cat(sprintf("%-7s", label),   sprintf("%7.3f", eval), "\n", sep="")
	cat(sprintf("%-7s", "Cont."), sprintf("%7.1f", eval/nv*100), "\n", sep="")
	cat(sprintf("%-7s", "Cum."),  sprintf("%7.1f", cumsum(eval/nv*100)), "\n", sep="")
}
# plot メソッド
plot.pfa <- function(	result,						# pfa が返すオブジェクト
			before=FALSE,					# 因子軸の回転前の結果を使うか回転後の結果を使うかを指定
			fac.no=1:2,					# 横軸と縦軸にプロットする因子番号
			scores=FALSE,					# 因子得点を描くか因子負荷量を描くかの指定
			xlab=NULL, ylab=NULL,				# 横軸，縦軸の名前
			axis=TRUE,					# 座標軸を描くかどうかの指定
			label.cex=0.7,					# 描画点につけるラベルの文字の大きさ
			...)						# plot 関数に引き渡す引数
{
	fac.name <- names(result$before.rotation.eval)
	if (length(fac.name) > 1) {
		ax1 <- fac.no[1]
		ax2 <- fac.no[2]
		if (is.null(xlab)) {
			xlab <- fac.name[ax1]
		}
		if (is.null(ylab)) {
			ylab <- fac.name[ax2]
		}
		if (scores) {
			x <- result$scores[, ax1]
			y <- result$scores[, ax2]
			labels <- 1:length(x)
		}
		else {
			if (before || is.na(result$after.rotation.fl)) {
				fl <- result$before.rotation.fl
			}
			else {
				fl <- result$after.rotation.fl
			}
			x <- fl[, ax1]
			y <- fl[, ax2]
			labels <- names(result$communality)
		}
		plot(x, y, xlab=xlab, ylab=ylab, ...)
		old <- par(xpd=TRUE)
		text(x, y, labels, cex=label.cex, pos=1)
		par(old)
	}
	if (axis) {
		abline(h=0, v=0)
	}
}
# φ係数を求める関数（chisq 関数が必要）
phi <- function(mat)	
{
	sqrt(chisq(mat)/sum(mat))
}

# コンティンジェンシー係数を求める関数（chisq 関数が必要）
contingency <- function(mat)	
{
	temp <- chisq(mat)
	sqrt(temp/(sum(mat)+temp))
}

# クラメール係数を求める関数（phi, chisq 関数が必要）
cramer <- function(mat)
{
	phi(mat)/sqrt(min(nrow(mat), ncol(mat))-1)
}

# カイ二乗値を計算する関数
chisq <- function(mat)
{
	ex <- outer(rowSums(mat), colSums(mat))/sum(mat)	# 期待値
	sum((mat-ex)^2/ex)					# カイ二乗値
}
# 図形描画のための関数群

# プロット関数群の開始。
# (x1, y1)-(x2, y2) 描画領域の左下隅と右上隅の座標を宣言する。
# 普通に R の関数（hist や plot などなど）を使って描かれたグラフに以下の関数を使って図形を付加する場合には不要。
plot.start <- function(x1=0, y1=0, x2=500, y2=500, ...)
{
	plot(c(x1, x2), c(y1, y2), type="n", xlab="", xaxt="n", ylab="", yaxt="n", bty="n", ...)
}

# (x1, y1)-(x2, y2) を結ぶ直線を描画する。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
plot.line <- function(x1, y1, x2, y2, ...)
{
	lines(c(x1, x2), c(y1, y2), ...)
}

# (x1, y)-(x2, y) を結ぶ水平な直線を描画する。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
plot.hline <- function(x1, x2, y, ...)
{
	lines(c(x1, x2), c(y, y), ...)
}

# (x, y1)-(x, y2) を結ぶ垂直な直線を描画する。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
plot.vline <- function(x, y1, y2, ...)
{
	lines(c(x, x), c(y1, y2), ...)
}

# (x1, y1)-(x2, y2) を対角頂点とする長方形（正方形）を描画する。
# 実際には polygon 関数を呼ぶので，... には polygon 関数が許容する引数を書くことができる。
plot.box <- function(x1, y1, x2, y2, ...)
{
	polygon(c(x1, x2, x2, x1), c(y1, y1, y2, y2), ...)
}

# 中心を (ox, oy) とする，半径 r の円を描く。
# start, end には，描き始めと描き終わりの位置を指定できる。
# 水平線を基準として，角度（度単位）で指定できる（0, 360 とすると，完全な円を描くことになる。
# 90, 270 とすると，左半分の円を描くことを指示することになる）。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
plot.circle <- function(ox, oy, r, start=0, end=360, ...)
{
	plot.ellipse(ox, oy, r, r, 0, start, end, ...)
}

# 中心を (ox, oy) とする，半径 r の円を描き，内部を塗りつぶす。
# 実際には polygon 関数を呼ぶので，... には polygon 関数が許容する引数を書くことができる。
plot.circlef <- function(ox, oy, r, ...)
{
	plot.ellipse(ox, oy, r, r, 0, 0, 360, func=polygon, ...)
}

# 度をラジアンに変換する
radian <- function(degree) {
	degree/180*pi
}

# 中心を (ox, oy) とする，長径 ra，短径 rb の楕円を描く。
# phi は楕円の傾き（長径が水平線となす角度），start, end には，描き始めと描き終わりの位置を指定できる。
# 長径を基準として，角度（度単位）で指定できる（0, 360 とすると，完全な楕円を描くことになる。
# 90, 270 とすると，左半分の楕円を描くことを指示することになる）。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
plot.ellipse <- function(ox, oy, ra, rb, phi=0, start=0, end=360, length=100, func=lines, ...)
{
	theta <- c(seq(radian(start), radian(end), length=length), radian(end))
	if (phi == 0) {
		func(ra*cos(theta)+ox, rb*sin(theta)+oy, ...)
	}
	else {
		x <- ra*cos(theta)
		y <- rb*sin(theta)
		phi <- radian(phi)
		cosine <- cos(phi)
		sine <- sin(phi)
		func(cosine*x-sine*y+ox, sine*x+cosine*y+oy, ...)
	}
}

# 中心を (ox, oy) とする，長径 ra，短径 rb の楕円を描き，内部を塗りつぶす。
# phi は楕円の傾き（長径が水平線となす角度）。
# 実際には polygon 関数を呼ぶので，... には polygon 関数が許容する引数を書くことができる。
plot.ellipsef <- function(ox, oy, ra, rb, phi=0, ...)
{
	plot.ellipse(ox, oy, ra, rb, phi, 0, 360, func=polygon, ...)
}

# (x1, y1) を描きはじめとして，右方向へ水平に l の長さの辺を描きその後反時計回りに辺を構成して n 正方角形を描く。
# 実際には polygon 関数を呼ぶので，... には polygon 関数が許容する引数を書くことができる。
plot.polygon <- function(x1, y1, l, n, ...)
{
	theta <- seq(0, 2*pi, length=n+1)
	x <- rep(x1, n)
	y <- rep(y1, n)
	for (i in 2:n) {
		x[i] <- x[i-1]+l*cos(theta[i])
		y[i] <- y[i-1]+l*sin(theta[i])
	}
	polygon(x, y, ...)
}

# (ox, oy) 中心とする円に内接する正 n 角形を描画する。
# 最初の頂点の位置は phi 引数（度）によって，反時計回りに再設定できる。
# 実際には polygon 関数を呼ぶので，... には polygon 関数が許容する引数を書くことができる。
plot.polygon2 <- function(ox, oy, r, n, phi=90, ...)
{
	theta <- seq(0, 2*pi, length=n+1)+radian(phi)
	polygon(r*cos(theta)+ox, r*sin(theta)+oy, ...)
}

# (x1, y1)-(x2, y2) を対角頂点とする長方形（正方形）を描画し，横方向の間隔が wx, 縦方向の間隔が wy となる格子を描く。
# 実際には lines 関数を呼ぶので，... には lines 関数が許容する引数を書くことができる。
# 長方形内部を塗りつぶすには，前もって plot.box で塗りつぶした長方形を描いてから plot.grid 関数を使う。
plot.grid <- function(x1, y1, x2, y2, wx, wy=NULL, ...)
{
	X1 <- min(x1, x2)
	X2 <- max(x1, x2)
	Y1 <- min(y1, y2)
	Y2 <- max(y1, y2)
	for (i in 0:as.integer(abs(X2-X1)/wx)) {
		plot.line(X1+wx*i, Y1, X1+wx*i, Y2, ...)
	}
	if (is.null(wy)) wy <- wx
	for (i in 0:as.integer(abs(y2-y1)/wy)) {
		plot.line(X1, Y1+wy*i, X2, Y1+wy*i, ...)
	}
}
# dual クラス のための plot メソッド　（dual, pc.dual, ro.dual が利用する）
plot.dual <- function(	x,					# dual が返すオブジェクト
			first=1,				# 横軸にプロットする解
			second=2,				# 縦軸にプロットする解
			weighted=FALSE,				# 相関比で重み付けした解をプロットするなら TRUE
			color.row="blue", color.col="black",	# 行と列のプロット色
			mark.row=19, mark.col=15,		# 行と列のプロット記号
			xlab=paste("Axis", first, sep="-"),	# 横座標軸名
			ylab=paste("Axis", second, sep="-"),	# 縦座標軸名
			axis=FALSE,				# 座標軸を点線で描くなら TRUE
			xcgx=FALSE,				# 横軸に取る座標の符号反転が必要なら TRUE 
			xcgy=FALSE,				# 縦軸に取る座標の符号反転が必要なら TRUE
			...)					# points, text 等に渡されるその他の引数
{
	if (ncol(x[[1]]) == 1) {
		warning("解が1個しかありません。二次元配置図は描けません。")
		return
	}
	suf <- if (weighted) 4 else 2				# 相関比で重み付けした解も選べる
	old <- par(xpd=TRUE, mar=c(5.1, 5.1, 2.1, 5.1))		# 左右を大きめに空ける
	row1 <- x[[suf]]  [, first]				# 横軸に取る解
	col1 <- x[[suf+1]][, first]
	if (xcgx) {						# 必要なら符号反転
		row1 <- -row1
		col1 <- -col1
	}
	row2 <- x[[suf]]  [, second]				# 縦軸に取る解
	col2 <- x[[suf+1]][, second]
	if (xcgy) {						# 必要なら符号反転
		row2 <- -row2
		col2 <- -col2
	}
	plot(c(row1, col1), c(row2, col2), type="n", xlab=xlab, ylab=ylab, ...)
	points(row1, row2, pch=mark.row, col=color.row, ...)
	text(row1, row2, labels=names(row1), pos=3, col=color.row, ...)
	points(col1, col2, pch=mark.col, col=color.col, ...)
	text(col1, col2, labels=names(col1), pos=3, col=color.col, ...)
	par(old)
	if (axis) {						# 座標軸を点線で描くならば
		abline(v=0, h=0, lty=3, ...)
	}
}
# ずっと以前には，R には plot.design 関数がなかったので，それに近いことをする関数を書いた
plot.design2 <- function(	data,				# データ行列
				FUN=mean)			# 作用させる関数
{
	nv <- ncol(data)					# データ行列の列数（変数の個数）
	result <-unlist(sapply(2:nv, function(i) by(data[,1], data[,i], FUN)))
	min.y <- min(result)
	max.y <- max(result)
	min.y <- min.y-(max.y-min.y)*0.1
	plot(c(0.5, nv-0.5), c(min.y, max.y), type="n", xlab="", xaxt="n", ylab=variable.names(data)[1])
	for (i in 2:nv) {
		r <- by(data[,1], data[,i], FUN)
		n <- length(r)
		lines(rep(i-1, n), r, type="o")
		text(i-1, min.y, variable.names(data)[i])
		text(i-1, r, paste("c", 1:n, sep=""), pos=4, offset=0.5)
	}
}
poisson.conf <- function(x, conf=0.95)
{
	N2 <- length(x)*2
	df <- 2*sum(x)+2
	alpha2 <- (1-conf)/2
	return(qchisq(c(alpha2, 1-alpha2), df)/N2)
}
# ポアソン分布への適合度の検定
poissondist <- function(d,                 # 度数ベクトル
                        x=NULL)            # 階級値ベクトル
{
  if (is.null(x)) {
    stop("関数の仕様が変更されました。度数ベクトルと同じ長さの階級値ベクトルも指定してください。")
  }
  data.name <- paste(deparse(substitute(d)), deparse(substitute(x)), sep=", ")
  method <- "ポアソン分布への適合度の検定"
  k <- length(d)                           # 階級数
  if (length(x) != k) {
    stop("度数ベクトル階級値ベクトルの長さは同じでなければなりません。")
  }
  	
  o <- numeric(diff(range(x))+1)
  o[x-min(x)+1] <- d
  x <- min(x):max(x)
  	
  k <- length(o)                           # 階級数
  n <- sum(o)                              # データ数
  lambda <- sum(o*x)/n                     # 平均値（=λ）

  p <- dpois(x, lambda)                    # 確率
  p[1] <- ppois(min(x), lambda)            # 最初と最後の階級値の確率は階級値以下・以上の確率を併合する
  p[k] <- ppois(max(x)-1, lambda, lower.tail=FALSE)
  e <- n*p                                 # 期待値
  table <- data.frame(x, o, p, e)          # 結果をデータフレームにまとめる
  rownames(table) <- paste("c-", x, sep="")

  while (e[1] < 1) {                       # 1 未満のカテゴリーの併合
    o[2] <- o[2]+o[1]
    e[2] <- e[2]+e[1]
    o <- o[-1]
    e <- e[-1]
    k <- k-1
  }
  while (e[k] < 1) {                       # 1 未満のカテゴリーの併合
    o[k-1] <- o[k-1]+o[k]
    e[k-1] <- e[k-1]+e[k]
    o <- o[-k]
    e <- e[-k]
    k <- k-1
  }
  chisq <- sum((o-e)^2/e)                  # カイ二乗統計量
  df <- k-2                                # 自由度
  p <- pchisq(chisq, df, lower.tail=FALSE) # P 値
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, parameter=df, p.value=p,
    estimate=c(n=n, lambda=lambda), method=method,
    data.name=data.name, table=table), class=c("htest", "poissondist")))
}
# summary メソッド
summary.poissondist <- function(obj,       # poissondist が返すオブジェクト
                                digits=5)
{
  table <- obj$table
  colnames(table) <- c("階級", "度数", "確率", "期待値")
  cat("\n適合度\n\n")
  print(table, digits=digits, row.names=FALSE)
}
# plot メソッド
plot.poissondist <- function(obj,          # poissondist が返すオブジェクト
                             ...)          # barplot へ渡す引数
{
  table <- obj$table
  nr <- nrow(table)
  posx <- barplot(table$o, space=0, ...)   # 観察度数を barplot で描く
  old <- par(xpd=TRUE)
  points(posx, table$e, pch=3)             # 理論度数を，記号 + で示す
  text(posx, -strheight("H"), table$x)
  par(old)
}
# 相関係数の検定・区間推定のパワーアナリシス
power.cor.test <- function(	n=NULL,						# 標本サイズ
				cor0=0,						# 母相関係数
				cor1=NULL,					# 標本相関係数
				sig.level=0.05,					# 有意水準
				power=NULL,					# 検出力
				alternative=c("two.sided", "one.sided"))	# 仮説・信頼区間の種類
{
	if (sum(sapply(list(n, cor0, cor1, power, sig.level), is.null)) != 1) {
		stop("en, delta, sd, power, sig.level のどれか一つだけを NULL にする")
	}
	alternative <- match.arg(alternative)					# 引数の補完
	tside <- switch(alternative, one.sided=1, two.sided=2)
	power.function <- quote(pnorm(sqrt(n-3)*abs(atanh(cor0)-atanh(cor1))-qnorm(sig.level/tside, lower.tail=FALSE)))
	if (is.null(power)) {
		power <- eval(power.function)
	}
	else if (is.null(n)) {
		n <- uniroot(function(n) eval(power.function)-power, c(4, 1e7))$root
	}
	else if (is.null(cor0)) {
		cor0 <- uniroot(function(cor0) eval(power.function)-power, c(cor1, 1))$root
	}
	else if (is.null(cor1)) {
		cor0 <- uniroot(function(cor1) eval(power.function)-power, c(cor0, 1))$root
	}
	else if (is.null(sig.level)) {
		sig.level <- uniroot(function(sig.level) eval(power.function)-power, c(1e-5, 0.99999))$root
	}
	else {
		stop("internal error")
	}
	METHOD <- "Power calculation of the one-sample correlation test."
	structure(list(n=n, cor0=cor0, cor1=cor1, sig.level=sig.level, power=power, alternative=alternative, method=METHOD), class="power.htest")
}
# 母比率の検定・区間推定のパワーアナリシス
power.one.sample.prop.test <- function(	n=NULL,						# 標本サイズ
					p=NULL,						# 標本比率
					p0=NULL,					# 母比率
					sig.level=0.05,					# 有意水準
					power=NULL,					# 検出力
					alternative=c("two.sided", "one.sided"))	# 仮説・信頼区間の種類
{
	if (sum(sapply(list(n, p, p0, power, sig.level), is.null)) != 1) {
		stop("n, p, p0, power, sig.level のどれか一つだけを NULL とする")
	}
	alternative <- match.arg(alternative)						# 引数の補完
	tside <- switch(alternative, one.sided=1, two.sided=2)
	power.function <- quote(pnorm((sqrt(n)*abs(p-p0)-qnorm(sig.level/tside, lower.tail=FALSE)*sqrt(p0*(1-p0)))/(sqrt(p*(1-p)))))
	if (is.null(power)) {
		power <- eval(power.function)
	}
	else if (is.null(n)) {
		n <- uniroot(function(n) eval(power.function)-power, c(1, 1e7))$root
	}
	else if (is.null(p)) {
		p <- uniroot(function(p) eval(power.function)-power, c(1e-5, p0))$root
	}
	else if (is.null(p0)) {
		p0 <- uniroot(function(p0) eval(power.function)-power, c(1e-5, p))$root
	}
	else if (is.null(sig.level)) {
		sig.level <- uniroot(function(sig.level) eval(power.function)-power, c(1e-5, 0.99999))$root
	}
	else {
		stop("internal error")
	}
	METHOD <- "Power calculation of the one-sample proportion test."
	structure(list(n=n, p=p, p0=p0, sig.level=sig.level, power=power, alternative=alternative, method=METHOD), class="power.htest")
}
# 母平均の検定・区間推定のパワーアナリシス
power.one.sample.t.test <- function(	n=NULL,					# 標本サイズ
					delta=NULL,				# 母平均との差
					sd=NULL,				# 母標準偏差
					sig.level=0.05,				# 有意水準
					power=NULL,				# 検出力
					alternative=c("two.sided", "one.sided"))# 仮説・信頼区間の種類
{
	if (sum(sapply(list(n, delta, sd, power, sig.level), is.null)) != 1) {
		stop("n, delta, sd, power, and sig.level のどれか一つだけが NULL でなければならない")
	}
	alternative <- match.arg(alternative)					# 省略された引数の補完
	tside <- switch(alternative, one.sided=1, two.sided=2)
	power.function <- quote(pnorm((sqrt(n)*delta/sd-qnorm(sig.level/tside, lower.tail=FALSE))))
	if (is.null(power)) {
		power <- eval(power.function)
	}
	else if (is.null(n)) {
		n <- uniroot(function(n) eval(power.function)-power, c(1, 1e7))$root
	}
	else if (is.null(delta)) {
		delta <- uniroot(function(delta) eval(power.function)-power, c(-1e7, 1e7))$root
	}
	else if (is.null(sd)) {
		sd <- uniroot(function(sd) eval(power.function)-power, c(1e-7, 1e7))$root
	}
	else if (is.null(sig.level)) {
		sig.level <- uniroot(function(sig.level) eval(power.function)-power, c(1e-5, 0.99999))$root
	}
	else {
		stop("internal error")
	}
	METHOD <- "1 標本 t 検定の検出力"
	structure(list(n=n, delta=delta, sd=sd, sig.level=sig.level, power=power, alternative=alternative, method=METHOD), class="power.htest")
}
# サンプルサイズが異なる二群の比率の差の検定の検出力
power.prop.test2 <- function(	Nc,					# 第一群のサンプルサイズ
				Nt,					# 第二群のサンプルサイズ
				Pc,					# 第一群での比率
				Pt,					# 第二群での比率
				sig.level=0.05)				# 有意水準
{
	P <- (Nc*Pc+Nt*Pt)/(Nc+Nt)					# プールした比率
	Z.alpha <- qnorm(sig.level/2, lower.tail=FALSE)
	Z.beta <- (abs(Pc-Pt)-Z.alpha*sqrt(P*(1-P)/Nc+P*(1-P)/Nt)) / sqrt(Pc*(1-Pc)/Nc+Pt*(1-Pt)/Nt)
	return(pnorm(Z.beta))						# 検出力
}
#=========================================================================================
# サンプルサイズが異なる二群の比率の差の検定の必要サンプルサイズ
power.prop.test3 <- function(	Pc,					# 第一群での比率
				Pt,					# 第二群での比率
				r=1,					# r = Nc/Nt
				sig.level=0.05,				# 有意水準
				power=0.8)				# 検出力
{
	P <- (r*Pc+Pt) / (r+1)						# プールした比率
	Z.alpha <- qnorm(sig.level/2, lower.tail=FALSE)
	Z.beta  <- qnorm(1-power, lower.tail=FALSE)
	Nt <- (Z.alpha*sqrt((r+1)*P*(1-P))+Z.beta*sqrt(r*Pt*(1-Pt)+Pc*(1-Pc)))^2 / r / (Pt-Pc)^2
	return(c(Nc=Nt*r, Nt=Nt))					# 第一群，第二群の必要サンプルサイズ
}
# 二群の平均値の差（両側検定）において，二群のサンプルサイズが異なるときの検出力を求める
power.t.test2 <- function(	n1,			# サンプルサイズ
				n2,			# サンプルサイズ
				delta,			# 効果量
				sig.level=0.05)		# 有意水準
{
	phi <- n1+n2-2					# 自由度
	lambda <- sqrt(n1*n2/(n1+n2))*delta
	q <- qt(sig.level/2, phi, lower.tail=FALSE)
	return(pt(-q, phi, ncp=lambda)+pt(q, phi, ncp=lambda, lower.tail=FALSE))
}
# prcomp を援用して，主成分分析の結果をまとめる
prcomp2 <- function(	dat,							# データ行列
			pcs=0,							# 求める主成分の個数
			cor=TRUE,						# TRUE の場合には相関係数行列，FALSE なら分散共分散行列を対象にする
			verbose=TRUE)						# 饒舌に結果を表示するかどうか
{
	p <- ncol(dat)								# 変数の個数
	if (is.null(colnames(dat))) colnames(dat) <- paste("Var", 1:p, sep=".")	# 変数名がついていないとき
	n <-nrow(dat)								# サンプルサイズ
	if (is.null(rownames(dat))) rownames(dat) <- paste("Obs", 1:n, sep=".")	# 行名がついていないとき
	dat <- subset(dat, complete.cases(dat))					# 欠損値を持つケースを除く	
	n <- nrow(dat)								# サンプルサイズは小さくなったかもしれない
	result<-prcomp(dat, scale=cor)						# prcomp を呼び出す
	if (pcs == 0) {								# 抽出する主成分数が指定されていないときは，
		pcs <- sum(result$sdev^2 >= 1)					# 1 以上の固有値の数
	}
	loadings <- t(result$sdev*t(result$rotation))[, 1:pcs, drop=FALSE]	# 主成分負荷量
	Contribution <- rowSums(loadings^2)					# 寄与率
	Eigen.values <- c(result$sdev[1:pcs]^2, sum(result$sdev[1:pcs]^2))	# 固有値
	denominator <- if(cor) p else sum(diag(cov(dat)))			# 累積寄与率を計算する分母（cor=FALSE なら，分散共分散行列の対角成分の和）
	Proportion <- Eigen.values/denominator*100				# 寄与率
	Cumulative.prop. <- cumsum(Proportion)					# 累積寄与率
	Cumulative.prop.[pcs+1] <- NA						# 末尾は不要
	result$loadings <- rbind(cbind(loadings, Contribution),			# 主要分析結果
					Eigen.values, Proportion, Cumulative.prop.)
	result$scores <- result$x * sqrt(n/(n-1))				# 主成分得点
	if (verbose) {
		cat(sprintf("%16s", ""), sprintf("%8s", paste("PC", 1:pcs, sep="")), " Contribution\n", sep="")
		for (i in 1:p) {
			cat(sprintf("%-16s", rownames(result$loadings)[i]), sprintf("%8.3f", loadings[i,]), sprintf("%10.3f\n", Contribution[i]), sep="")
		}
		cat("Eigen.values    ", sprintf("%8.3f", result$loadings[p+1, 1:pcs]), "\n", sep="")
		cat("Proportion      ", sprintf("%8.3f", result$loadings[p+2, 1:pcs]), "\n", sep="")
		cat("Cumulative.prop.", sprintf("%8.3f", result$loadings[p+3, 1:pcs]), "\n", sep="")
	}
	invisible(result)
}
# 主座標分析を行う
princo <- function(	s)				# 類似度行列
{
	if (is.data.frame(s)) {
		s <- as.matrix(s)
	}
	n <- nrow(s)					# 行列のサイズ
	object.names <- colnames(s)
	if (is.null(object.names)) {
		object.names <- paste("対象", 1:n, sep="")
	}
	m <- colSums(s)/n
	h <- s+sum(s)/n^2-outer(m, m, "+")
	res <- eigen(h)					# 固有値・固有ベクトルを求める
	values <- res$values[res$values > 1e-5]		# 解の個数を決める
	ax <- length(values)
	vectors <- t(t(res$vectors[,1:ax])*sqrt(values))# 対象に数値を与える
	colnames(vectors) <- names(values) <- paste("解", 1:ax, sep="")
	rownames(vectors) <- object.names
	return(structure(list(ax=ax, n=n, values=values, vectors=vectors), class="princo"))
}
# print メソッド
print.princo <- function(	res,			# princo が返すオブジェクト
				ax=res$ax,		# 何次元までの解を出力するか
				digits=5)		# 表示桁数
{
	ax <- min(ax, res$ax)
	val <- res$values
	val2 <- val/sum(val)
	val <- rbind(val, val2, cumsum(val2))
	rownames(val) <- c("固有値", "寄与率", "累積寄与率")
	print(val[, 1:ax], digits=digits)
	cat("\nベクトル\n\n")
	print(res$vectors[, 1:ax], digits=digits)
}
# plot メソッド
plot.princo <- function(res,				# princo が返すオブジェクト
			labels=FALSE,			# ラベルを描くかどうか
			text.cex=0.7,			# ラベルのフォントの大きさ
			...)				# plot への引数
{
	if (res$ax >= 2) {				# 二次元以上の解が得られたら，
		plot(res$vectors[,1:2], ...)		# 二次元の図を描く
		abline(v=0, h=0)
		old <- par(xpd=TRUE)
		if (labels) {
			text(res$vectors[,1], res$vectors[,2], rownames(res$vectors),
			     pos=4, offset=.2, cex=text.cex)
		}
		par(old)
	}
	else {
		warning("解が一次元なので，二次元配置図は描けません。")
	}
}
# princomp を援用して，主成分分析の結果をまとめる
princomp2 <- function(	dat,							# データ行列
			pcs=0,							# 求める主成分の個数
			cor=TRUE,						# TRUE の場合には相関係数行列，FALSE なら分散共分散行列を対象にする
			scores=FALSE,						# 主成分得点を計算するかどうか
			verbose=TRUE)						# 饒舌に結果を表示するかどうか
{
	p <- ncol(dat)								# 変数の個数
	if (is.null(colnames(dat))) colnames(dat) <- paste("Var", 1:p, sep=".")	# 変数名がついていないとき
	n <- nrow(dat)								# サンプルサイズ
	if (is.null(rownames(dat))) rownames(dat) <- paste("Obs", 1:n, sep=".")	# 行名がついていないとき
	dat <- subset(dat, complete.cases(dat))					# 欠損値を持つケースを除く
	n <- nrow(dat)								# サンプルサイズは小さくなったかもしれない
	result<-princomp(dat, cor=cor, scores=scores)				# princomp を呼び出す
	if (pcs == 0) {								# 抽出する主成分数が指定されていないときは，
		pcs <- sum(result$sdev^2 >= 1)					# 1 以上の固有値の数
	}
	loadings <- t(result$sde*t(result$loadings))[, 1:pcs, drop=FALSE]	# 主成分負荷量
	Contribution <- rowSums(loadings^2)					# 寄与率
	Eigen.values <- c(result$sdev[1:pcs]^2, sum(result$sdev[1:pcs]^2))	# 固有値
	denominator <- if (cor) p else sum(diag(var(dat)*(n-1)/n))		# 累積寄与率を計算する分母（cor=FALSE なら，分散共分散行列の対角成分の和）
	Proportion <- Eigen.values/denominator*100				# 寄与率
	Cumulative.prop. <- cumsum(Proportion)					# 累積寄与率
	Cumulative.prop.[pcs+1] <- NA						# 末尾は不要
	result$loadings <- rbind(cbind(loadings, Contribution),			# 主要分析結果
				Eigen.values, Proportion, Cumulative.prop.)
	result$x <- result$scores * sqrt((n-1)/n)				# 主成分得点
	if (verbose) {
		cat(sprintf("%16s", ""), sprintf("%8s", paste("PC", 1:pcs, sep="")), " Contribution\n", sep="")
		for (i in 1:p) {
			cat(sprintf("%-16s", rownames(result$loadings)[i]), sprintf("%8.3f", loadings[i,]), sprintf("%10.3f\n", Contribution[i]), sep="")
		}
		cat("Eigen.values    ", sprintf("%8.3f", result$loadings[p+1, 1:pcs]), "\n", sep="")
		cat("Proportion      ", sprintf("%8.3f", result$loadings[p+2, 1:pcs]), "\n", sep="")
		cat("Cumulative.prop.", sprintf("%8.3f", result$loadings[p+3, 1:pcs]), "\n", sep="")
		if (scores) {
			cat(sprintf("\n%16s", ""), sprintf("%8s", paste("PC", 1:pcs, sep="")), "\n", sep="")
			for (i in 1:n) {
				cat(sprintf("%-16s", rownames(result$scores)[i]), sprintf("%8.3f", result$scores[i, 1:pcs]), "\n", sep="")
			}
		}
	}
	invisible(result)
}
# ANOVA 表のプリント・メソッド（SAB, ASB 関数が返すオブジェクトを表示する）
print.anova.table <- function(x)
{
	printf <- function(x, fmt) if (is.na(x)) "" else sprintf(fmt, x)
	x[,4] <- sapply(x[,4], printf, "%.5f")
	x[,5] <- sapply(x[,5], printf, "%.5f")
	print.data.frame(x)
}
# データフレーム・行列をタブ区切りで整形して書き出す
print.fixed <- function(x,				# データフレームまたは行列。クラスとして "fixed" を与えるとよいこともある
			format=NULL)			# 表示書式。format="f5 2i s" など
{
	if (!is.matrix(x)) {				# 行列でないなら，
		x <- as.matrix(x)			# 行列にしてしまう
	}
	nr <- nrow(x)					# 行数
	nc <- ncol(x)					# 列数
	if (is.null(format)) {				# format が NULL なら，
		format <- rep("%20s", nc)		# 適当に作る
	}
	else {						# format が与えられているなら，
		format <- unlist(strsplit(format, "[ ,]"))	# sprintf で使える，正式な format に変換する
		format <- unlist(sapply(format, function(str) {	# Rf3 のような形式（R 繰り返し）に対応
			rep <- sub("[a-z]+[0-9]*", "", str)
			if (rep == "") {
				rep <- "1"
			}
			rep <- as.integer(rep)
			fmt <- sub("^[0-9]*", "", str)
			return(rep(fmt, rep))
		}))
		if (nc != length(format)) {		# 書式は十分か？
			stop("insufficient format")
		}
		format <- sapply(1:nc, function(i) {
			fmt <- format[i]
			ch <- substring(fmt, 1, 1)
			if (ch == "f") {
				format[i] <- paste("%.", substring(fmt, 2), "f", sep="")
			}
			else {
				format[i] <- paste("%", ch, sep="")
			}
		})
	}
	for (i in 1:nr) {
		for (j in 1:nc) {
			number <- suppressWarnings(as.numeric(x[i,j]))
			if (is.na(number)) {
				if (x[i,j] == "NA") x[i,j] = ""
				cat(sprintf("%s", x[i,j]), if (j==nc) "\n" else "\t", sep="")
			}
			else {
				cat(sprintf(format[j], number), if (j==nc) "\n" else "\t", sep="")
			}
		}
	}
}
print.htest <- function (x, digits = 4, quote = TRUE, prefix = "", ...) 
{
    conv1 <- function(x) ##### 検定の名称
    {
#     	print(x)
	if(any(grep("[0-9]+-sample test for equality of proportions without continuity correction", x))) {
		x <- "割合の一様性の検定（連続性の補正なし）"
	}
    	else if(any(grep("[0-9]+-sample proportions test without continuity correction", x))) {
		x <- "割合の検定（連続性の補正なし）"
	}
    	else if (x == "Box-Pierce test") { # Box.test
    		x <- "Box-Pierce検定"
    	}
    	else if (x == "Box-Ljung test") {
    		x <- "Box-Ljung検定"
    	}
    	else if (x == "Phillips-Perron Unit Root Test") { # PP.test
    		x <- "Phillips-Perronの単位根検定"
    	}
    	else if (x == "Ansari-Bradley test") { # ansari.test
    		x <- "アンサリ・ブラドレイ検定"
    	}
     	else if (x == "Bartlett test of homogeneity of variances") { # bartlett.tes
    		x <- "分散の一様性の検定（バートレット検定）"
    	}
   	else if (x == "Exact binomial test") { # binom.test
    		x <- "二項検定"
    	}
    	else if (x == "Pearson's Chi-squared test with Yates' continuity correction") { # chisq.test
    		x <- "ピアソンのカイ二乗検定（イエーツの連続性補正）"
    	}
     	else if (x == "Pearson's Chi-squared test") { # chisq.test
    		x <- "ピアソンのカイ二乗検定（連続性補正なし）"
    	}
   	else if (x == "Chi-squared test for given probabilities") { # chisq.test
    		x <- "理論比が与えられたときのカイ二乗検定（適合度検定）"
    	}
    	else if (x == "Pearson's product-moment correlation") { # cor.test
    		x <- "ピアソンの積率相関係数"
    	}
    	else if (x == "Spearman's rank correlation rho") { # cor.test
    		x <- "スピアマンの順位相関係数"
    	}
    	else if (x == "Kendall's rank correlation tau") { # cor.test
    		x <- "ケンドールの順位相関係数"
    	}
    	else if (x == "Fisher's Exact Test for Count Data") { # fisher.test
    		x <- "計数データにおけるフィッシャーの正確確率検定"
    	}
    	else if (x == "Fligner-Killeen test of homogeneity of variances") { # fligner.test
    		x <- "分散の一様性の検定（Fligner-Killeen検定）"
    	}
    	else if (x == "Friedman rank sum test") { # fredman.test
    		x <- "フリードマン検定"
    	}
    	else if (x == "Kruskal-Wallis rank sum test") { # kruskal.test
    		x <- "クラスカル・ウォリス検定"
    	}
    	else if (x == "Two-sample Kolmogorov-Smirnov test") { # ks.test
    		x <- "二標本コルモゴロフ・スミルノフ検定"
    	}
    	else if (x == "One-sample Kolmogorov-Smirnov test") { # ks.test
    		x <- "一標本コルモゴロフ・スミルノフ検定"
    	}
    	else if (x == "Mantel-Haenszel chi-squared test with continuity correction") { # mantelhaen.test
    		x <- "マンテル・ヘンツェルのカイ二乗検定（連続性の補正）"
    	}
    	else if (x == "Exact conditional test of independence in 2 x 2 x k tables") { # mantelhaen.test
    		x <- "2 x 2 x k 分割表における条件付き独立性の正確検定"
    	}
    	else if (x == "Cochran-Mantel-Haenszel test") { # mantelhaen.test
    		x <- "コクラン・マンテル・ヘンツェル検定"
    	}
    	else if (x == "McNemar's Chi-squared test with continuity correction") { # mcnemar.test
    		x <- "マクネマー検定（連続性の補正）"
    	}
    	else if (x == "McNemar's Chi-squared test") { # mcnemar.test
    		x <- "マクネマー検定（連続性の補正なし）"
    	}
    	else if (x == "Mood two-sample test of scale") { # mood.test
    		x <- "尺度についての二標本ムード検定"
    	}
     	else if (x == "One-way analysis of means") { # oneway.test
    		x <- "一元配置分散分析"
    	}
     	else if (x == "One-way analysis of means (not assuming equal variances)") { # oneway.test
    		x <- "一元配置分散分析（等分散を仮定しない場合）"
    	}
     	else if (x == "Chi-squared Test for Trend in Proportions") { # prop.trend.test
    		x <- "割合の傾向についてのカイ二乗検定（傾向検定）"
    	}
     	else if (x == "Quade test") { # quade.test
    		x <- "Quade検定"
    	}
     	else if (x == "Shapiro-Wilk normality test") { # shapiro.test
    		x <- "シャピロ・ウィルクの正規性検定"
    	}
    	else if (x == "Welch Two Sample t-test") { # t.test
    		x <- "二標本t検定（Welchの方法）"
    	}
    	else if (x ==" Two Sample t-test") { # t.test
    		x <- "二標本t検定（分散が等しいと仮定できるとき）"
    	}
    	else if (x == "One Sample t-test") { # t.test
    		x <- "一標本t検定（母平均の検定）"
    	}
    	else if (x == "Paired t-test") { # t.test
    		x <- "対応のある場合のt検定"
    	}
    	else if (x == "F test to compare two variances") { # var.test
    		x <- "二群の等分散性の検定"
    	}
    	else if (x == "Wilcoxon signed rank test") { # wilcox.test
    		x <- "ウィルコクソンの符号付順位和検定"
    	}
    	else if (x == "Wilcoxon rank sum test") { # wilcox.test
    		x <- "ウィルコクソンの順位和検定（マン・ホイットニーのU検定）"
    	}
    	else if (x == "Wilcoxon rank sum test with continuity correction") { # wilcox.test
    		x <- "ウィルコクソンの順位和検定（連続性の補正）"
    	}
    	else if (x == "Wilcoxon signed rank test with continuity correction") { # wilcox.test
    		x <- "ウィルコクソンの符号付順位和検定（連続性の補正）"
    	}
    	return(x)
    }
    conv2 <- function(x)
    {
    	if (length(x) == 2) {
    		if (x[1] == "num df") {
    			x[1] <- "第1自由度"
    		}
    		if (x[2] == "denom df") {
    			x[2] <- "第2自由度"
    		}
    	}
    	else if (x == "df") {
    		x <- "自由度"
    	}
    	else if (x == "Truncation lag parameter") {
    		x <- "切り捨てラグ・パラメータ"
    	}
    	else if (x == "number of trials") {
    		x <- "試行数"
    	}
    	return(x)
    }
    conv3 <- function(x) ##### 検定対象の名前
    {
    	if (x =="difference in means") {
    		x <- "母平均の差"
    	}
    	else if (x == "mean") {
    		x <- "母平均"
    	}
    	else if (x == "mu") {
    		x <- "母平均"
    	}
    	else if (x == "correlation") {
    		x <- "母相関"
    	}
    	else if (x == "rho") {
    		x <- "母相関（ρ）"
    	}
    	else if (x == "tau") {
    		x <- "母相関（τ）"
    	}
    	else if (x == "probability of success") {
    		x <- "成功確率（母比率）"
    	}
    	else if (x == "ratio of scales") {
    		x <- "尺度の比"
    	}
    	else if (x == "common odds ratio") {
    		x <- "共通オッズ比"
    	}
    	else if (x == "odds ratio") {
    		x <- "オッズ比"
    	}
    	else if (x == "ratio of variances") {
    		x <- "分散比"
    	}
    	return(x)
    }
    conv4 <- function(x) ##### 推定値の名前
    {
    	names(x) <- sub("mean of the differences", "差の平均値", names(x))
    	names(x) <- sub("mean of ", "平均値", names(x))
    	names(x) <- sub("prop ", "割合", names(x))
     	if (any(grep("mean in group [0-9]+", names(x)[1]))) {
    		names(x) <- paste("グループ", 1:length(x), "の平均値", sep="")
    	}
    	names(x) <- sub("cor", "相関係数", names(x))
    	names(x) <- sub("rho", "ρ", names(x))
    	names(x) <- sub("tau", "τ", names(x))
    	names(x) <- sub("probability of success", "成功確率（母比率）", names(x))
    	names(x) <- sub("ratio of scales", "尺度の比", names(x))
    	names(x) <- sub("common odds ratio", "共通オッズ比", names(x))
    	names(x) <- sub("odds ratio", "オッズ比", names(x))
    	names(x) <- sub("ratio of variances", "分散比", names(x))
    	names(x) <- sub("difference in location", "位置母数の差", names(x))
    }
    conv5 <- function(x)
    {
    	if (any(grep(" and ", x))) {
    		return(gsub(" and ", " と ", x))
    	}
    	else if (any(grep(" by ", x))) {
    		y <- unlist(strsplit(x, " "))
    		return(paste(y[1], "を", y[3], "で層別"))
    	}
    	else if (any(grep("null probability", x))) {
    		return(gsub("null probability", "帰無仮説における母比率", x))
    	}
    	else if (any(grep("using scores", x))) {
    		return(gsub("using scores", "使用したスコア", x))
    	}
    	else {
    		return(x)
    	}
    }
    conv6 <- function(x) ##### 検定統計量の名前 STATISTIC
    {
    	if (length(x) == 1) {
	    	if (x == "X-squared") {
	    		x <- "カイ二乗値"
	    	}
	    	else if (x == "t") {
	    		x <- "t値"
	    	}
	    	else if (x == "Z") {
	    		x <- "Z値"
	    	}
	    	else if (x == "Bartlett's K-squared") {
	    		x <- "バートレットのK二乗値"
	    	}
	    	else if (x == "number of successes") {
	    		x <- "成功数"
	    	}
	    	else if (x == "Friedman chi-squared") {
	    		x <- "フリードマンのカイ二乗値"
	    	}
	    	else if (x == "Mantel-Haenszel X-squared") {
	    		x <- "マンテル・ヘンツェルのカイ二乗値"
	    	}
	    	else if (x == "Cochran-Mantel-Haenszel M^2") {
	    		x <- "コクラン・マンテル・ヘンツェルのM二乗値"
	    	}
	    	else if (x == "Fligner-Killeen:med chi-squared") {
	    		x <- "Fligner-Killeenのカイ二乗値"
	    	}
	    	else if (x == "Kruskal-Wallis chi-squared") {
	    		x <- "クラスカル・ウォリスのカイ二乗値"
	    	}
	    	else if (x == "McNemar's chi-squared") {
	    		x <- "マクネマーのカイ二乗値"
	    	}
    	}
    	return(x)
    }
  #  cat("\n")
    #@1
    writeLines(strwrap(conv1(x$method)))
    cat("\n")
    #@5
    cat("データ: ", conv5(x$data.name), "\n")
    out <- character()
    if (!is.null(x$statistic))
        #@6
        out <- c(out, paste(conv6(names(x$statistic)), "=", format(round(x$statistic, 
            4))))
    if (!is.null(x$parameter))
        #@2
        out <- c(out, paste(conv2(names(x$parameter)), "=", format(round(x$parameter, 
            3))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        out <- c(out, paste("P値", if (substr(fp, 1, 1) == 
            "<") fp else paste("=", fp)))
    }
    writeLines(strwrap(paste(out, collapse = ", ")))
    if (!is.null(x$alternative)) {
        cat("対立仮説: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1) {
                alt.char <- switch(x$alternative, two.sided = "ではない", 
                  less = "より小さい", greater = "より大きい")
                #@3
                cat(conv3(names(x$null.value)), "は，", x$null.value, alt.char, 
                  "\n", sep="")
            }
            else {
                cat(x$alternative, "\nnull values:\n")
                print(x$null.value, ...)
            }
        }
        else {
        	alt.char <- switch(x$alternative, two.sided = "等しくない", 
                  less = "小さい", greater = "大きい")
        	cat(alt.char, "\n")
        }
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), "パーセント信頼区間: ", 
            format(c(x$conf.int[1], x$conf.int[2])), "\n")
    }
    if (!is.null(x$estimate)) {
        cat("標本推定値: \n")
        #@4
        names(x$estimate) <- conv4(x$estimate)
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}
# データフレーム・行列を LaTeX の表にする
print.latex <- function(
			x, 			# matrix or data frame;  "latex" class
			file="", 		# output file (default to screen)
			append=FALSE, 		# rewrite or append output file
			ctable=TRUE, 		# use ctable.sty or not
			caption="", 		# caption
			cap="", 		# short caption (ctable only)
			label="", 		# label
			htbp="htbp", 		# layout			ex. htpb="ht" or "p" etc.
			align=NULL, 		# algnment of column		ex. align="rlcrlcrlc"
			format=NULL,		# format of column		ex. format="f5 2i s" == "%.5f %i %i %s"
			tnote=NULL, 		# tnote (ctable only)		ex. tnote=c("xxx...", "yyy...",...)
			line.feed=NULL, 	# \\[nn] (not for ctable)	ex. line.feed="-1.5mm" etc. 
			na2blank=TRUE		# display blank for NA values
		)
{
	nc <- ncol(x)									# オブジェクトの列数
	nr <- nrow(x)									# オブジェクトの行数
	has.rowname <- !is.null(rownames(x))						# 行名がある
	if (is.null(colnames(x))) {							# 列名がないときは，
		colnames(x) <- paste("Col", 1:nc, sep="-")				# "Col-##" の形式で名前を付ける
	}
	if (is.null(align)) {								# align の指定がないときは，
		align <- paste(rep("c", nc+has.rowname), collapse="")			# センター揃え "c…" にする
	}
	if (is.null(format)) {								# format の指定がないときは，
		format <- rep("%s", nc)							# "%s" で書く
	}
	else {										# format の指定があるときは，
		format <- unlist(strsplit(format, "[ ,]"))				# 文字ベクトルに分解
		format <- unlist(sapply(format, function(str) {				# Rf3 のような形式（R 繰り返し）に対応
			rep <- sub("[a-z]+[0-9]*", "", str)
			if (rep == "") {
				rep <- "1"
			}
			rep <- as.integer(rep)
			fmt <- sub("^[0-9]*", "", str)
			return(rep(fmt, rep))
		}))
		if (nc+has.rowname > length(format)) stop("too few formats")		# 行名がないときは列数に等しいこと，行名があれば列数プラス1
		else if (nc+has.rowname < length(format)) stop("too many formats")
		format <- sapply(1:(nc+has.rowname), function(i) {			# 正式な format に変換する
			fmt <- format[i]
			ch <- substring(fmt, 1, 1)
			if (ch == "f") {
				format[i] <- paste("%.", substring(fmt, 2), "f", sep="")
			}
			else {
				format[i] <- paste("%", ch, sep="")
			}
		})
	}
	cat("%\n", file=file, append=append, sep="")					# ヘッダー
	if (ctable) {									# ctable 形式で出力する場合
		cat(	"\\ctable[\n",
			"\tcap = {", cap, "},\n",
			"\tcaption = {", caption, "},\n",
			"\tlabel = ", label, ",\n",
			"\tpos = htbp\n",
			"]\n",
			"{", align, "}\n",
			file=file, append=append, sep="")
		if (is.null(tnote)) {
			cat("{}\n", file=file, append=append, sep="")
		}
		else {
			cat("{\n", file=file, append=append, sep="")
			sapply(1:length(tnote),
				function(i)
				cat("\\tnote[", letters[i], "]{", tnote[i], "}\n", file=file, append=append, sep=""))
			cat("}\n", file=file, append=append, sep="")
		}
		cat("{ \\FL\n", file=file, append=append, sep="")
		head <- sprintf("\t %s %s \\ML\n", if (has.rowname) "&" else "", paste(colnames(x), collapse=" & "))
		cat(head, file=file, append=TRUE)
		for (i in 1:nr) {
			temp <- sapply(1:nc, function(j) sprintf(format[j+has.rowname], x[i, j]))
			if (has.rowname) {
				body <- sprintf("\t %s & %s", sprintf(format[1], rownames(x)[i]), paste(temp, collapse=" & "))
			}
			else {
				body <- sprintf("\t %s", paste(temp, collapse=" & "))
			}
			if (na2blank) {							# na2blank == TRUE なら，
				body <- gsub("(^NA | NA | NA$)", " ", body)		# NA を出力しない（空白を出力する）
			}
			if (i == nr) {
				body <- sprintf("%s \\LL\n", body)
			}
			else {
				body <- sprintf("%s \\NN\n", body)
			}
			cat(body, file=file, append=TRUE)
		}
		cat("}\n", file=file, append=TRUE, sep="")
	}
	else {										# table 形式で出力する場合
		cat(	"\\begin{table}[", htbp, "]\n",
			"\t\\caption{", caption, "}\n",
			"\t\\label{", label, "}\n",
			"\t\\centering\n",
			"\t\t\\begin{tabular}{", align, "} \\hline\n",
			file=file, append=append, sep="")
		head <- sprintf("\t\t\t\t %s %s \\\\ \\hline\n", if (has.rowname) "&" else "", paste(colnames(x), collapse=" & "))
		cat(head, file=file, append=TRUE)
		for (i in 1:nr) {
			temp <- sapply(1:nc, function(j) sprintf(format[j+has.rowname], x[i, j]))
			if (has.rowname) {
				body <- sprintf("\t\t\t %s & %s",  sprintf(format[1], rownames(x)[i]),  paste(temp, collapse=" & "))
			}
			else {
				body <- sprintf("\t\t\t %s", paste(temp, collapse=" & "))
			}
			if (na2blank) {							# na2blank == TRUE なら，
				body <- gsub("(^NA | NA | NA$)", " ", body)		# NA を出力しない（空白を出力する）
			}
			if (i == nr) {
				body <- sprintf("%s \\\\ \\hline\n", body)
			}
			else if (is.null(line.feed)) {
				body <- sprintf("%s \\\\\n", body)
			}
			else {
				body <- sprintf("%s \\\\[%s]\n", body, line.feed)
			}
			cat(body, file=file, append=TRUE)
		}
		cat("\t\t\\end{tabular}\n", "\\end{table}\n", file=file, append=TRUE, sep="")
	}
	cat("%\n", file=file, append=TRUE, sep="")					# フッター
}
# C 言語の printf をシミュレートする
printf <- function(fmt, ...)
{
	cat(sprintf(fmt, ...))
}
# 母比率の信頼区間
prop.conf <- function(	r,				# 標本のうち，注目している特性を持つものの数
			n,				# サンプルサイズ
			conf=0.95,			# 信頼率（信頼度）
			approximation=FALSE)		# 正規分布による近似を行う場合に TRUE を指定する
{
	p <- r/n					# 標本比率
	alpha <- 1-conf					# α
	if (p == 0) {					# 標本比率が 0 の場合
		pl <- 0
		pu <- 1-alpha^(1/n)
	}
	else if (p == 1) {				# 標本比率が 1 の場合
		pl <- alpha^(1/n)
		pu <- 1
	}
	else if (approximation) {			# 正規分布による近似を行う場合
		z <- qnorm(alpha/2, lower.tail=FALSE)	# 両側確率がαになる Z 値
		x <- n/(n+z^2)*(p+z^2/(2*n)+c(1, -1)*z*sqrt(p*(1-p)/n+z^2/(4*n^2)))
		pu <- x[1]
		pl <- x[2]
	}
	else {						# F 分布を使って正確な信頼区間を求める場合
		nu1 <- 2*(n-r+1)
		nu2 <- 2*r
		Fv <- qf(alpha/2, nu1, nu2, lower.tail=FALSE)
		pl <- nu2/(nu1*Fv+nu2)
		nu1 <- 2*(r+1)
		nu2 <- 2*(n-r)
		Fv <- qf(alpha/2, nu1, nu2, lower.tail=FALSE)
		pu <- nu1*Fv/(nu1*Fv+nu2)
	}
	setNames(c(pl, pu), c("lower.p", "upper.p"))
}
# ものぐさ太郎，御用達の関数
qq <- function() q("no")
# 数量化 I 類
qt1 <- function(dat,						# データ行列
		y,						# 従属変数
		func.name=c("solve", "ginv"))			# 逆行列を取る関数名の選択
{
	vname <- colnames(dat)					# 変数名
	vname.y <- deparse(substitute(y))
	cname <- unlist(sapply(dat, levels))			# カテゴリー名
	dat <- data.frame(dat,  y)
	dat <- subset(dat, complete.cases(dat))			# 欠損値を持つケースを除く
	p <- ncol(dat)						# 最右端が従属変数，残りはカテゴリー変数
	ncat <- p-1						# 最右端が従属変数，残りはカテゴリー変数
	stopifnot(all(sapply(dat[, 1:ncat], is.factor)))	# 独立変数はすべて factor であること
	dat[, 1:ncat] <- lapply(dat[ ,1:ncat, drop=FALSE], as.integer)
	nc <- nrow(dat)						# 列数
	mx <- sapply(dat[, 1:ncat, drop=FALSE], max)		# 各カテゴリー変数の取る値の最大値
	start <- c(0, cumsum(mx)[-ncat])			# 延べの順序番号
	nobe <- sum(mx)						# 延べカテゴリー数

	# ダミー変数を使ったデータ行列に変換
	x <- t(apply(dat, 1,
			function(obs)
			{
				zeros <- numeric(nobe)
				zeros[start+obs[1:ncat]] <- 1
				c(zeros[-start-1], obs[ncat+1])
			}
		))

	# 重回帰分析として解く	
	a <- cov(x)
	ndim <- nobe-ncat
	if (match.arg(func.name) == "solve") {
		inverse <- solve
		B <- inverse(a[1:ndim, 1:ndim], a[ndim+1, 1:ndim])
	}
	else {
		library(MASS)
		inverse <- ginv
		B <- inverse(a[1:ndim, 1:ndim]) %*% a[ndim+1, 1:ndim]
	}
	m <- colMeans(x)
	const <- m[ndim+1]-sum(B*m[1:ndim])
	prediction <- x[,1:ndim]%*%as.matrix(B)+const
	observed <- x[,ndim+1]
	prediction <- cbind(observed, prediction, observed-prediction)

	# 数量化 I 類としての解に変換
	ncase <- nrow(dat)
	s <- colSums(x)
	name <- coef <- NULL
	en <- 0
	for (i in 1:ncat) {
		st <- en+1
		en <- st+mx[i]-2
		target <- st:en
		temp.mean <- sum(s[target]*B[target])/ncase
		const <- const+temp.mean
		coef <- c(coef, -temp.mean, B[target]-temp.mean)
	}
	coef <- c(coef, const)
	names(coef) <- c(paste(rep(vname, mx), cname, sep="."), "定数項")

	# アイテム変数と従属変数間の偏相関係数
	par <- matrix(0, nrow=nc, ncol=ncat)
	for (j in 1:nc) {
		en <- 0
		for (i in 1:ncat) {
			st <- en+1
			en <- st+mx[i]-2
			target <- st:en
			par[j, i] <- crossprod(x[j, target], B[target])
		}
	}
	par <- cbind(par, observed)
	r <- cor(par)
	print(vname)
	i <- inverse(r)
	d <- diag(i)
	partial.cor <- (-i/sqrt(outer(d, d)))[ncat+1, 1:ncat]
	partial.t <- abs(partial.cor)*sqrt((nc-ncat-1)/(1-partial.cor^2))
	partial.p <- pt(partial.t, nc-ncat-1, lower.tail=FALSE)*2
	partial <- cbind(partial.cor, partial.t, partial.p)

	coef <- as.matrix(coef)
	colnames(coef) <- "カテゴリースコア"
	colnames(prediction) <- c("観察値", "予測値", "残差")
	colnames(partial) <- c("偏相関係数", "t 値", "P 値")
	rownames(prediction) <- paste("#", 1:nc, sep="")
	rownames(partial) <- vname
	rownames(r) <- colnames(r) <- c(vname, vname.y)
	return(structure(list(coefficients=as.matrix(coef),
		r=r, partial=partial, prediction=prediction), class="qt1"))
}
# print メソッド
print.qt1 <- function(	obj,					# qt1 が返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	print(round(obj$coefficients, digits=digits))
}
# summary メソッド
summary.qt1 <- function(obj,					# qt1 が返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.qt1 <- function(	obj,					# qt1 が返すオブジェクト
			which=c("category.score", "fitness"),	# カテゴリースコアを表示するか観察値と予測値を表示するかの選択
			...)					# barplot, plot に引き渡す引数
{
	if (match.arg(which) == "category.score") {
		coefficients <- obj$coefficients[-length(obj$coefficients),]
		coefficients <- rev(coefficients)
		cname <- names(coefficients)
		names(coefficients) <- NULL
		barplot(coefficients, horiz=TRUE, xlab="カテゴリースコア", ...)
		text(0, 1.2*(1:length(cname)-0.5), cname, pos=ifelse(coefficients > 0, 2, 4))
	}
	else {
		result <- obj$prediction
		plot(result[, 2], result[, 1],  xlab="予測値", ylab="観察値", asp=1, ...)
		abline(c(0,1))
	}
}
# 数量化 II 類
qt2 <- function(dat,						# データ行列（アイテムデータ）
	group)		# 群変数
{
	geneig <- function(a, b)				# 一般化固有値問題を解く関数
	{
	    a <- as.matrix(a)
	    b <- as.matrix(b)
	    if (nrow(a) == 1) {
	        res <- list(values=as.vector(a/b), vectors=as.matrix(1))
	    }
	    else {
	        res <- eigen(b)
	        g <- diag(1/sqrt(res$values))
	        v <- res$vectors
	        res <- eigen(g %*% t(v) %*% a %*% v %*% g)
	        res$vectors <-v %*% g %*% res$vectors
	    }
	    return(res)
	}

	ok <- complete.cases(dat, group)
	dat <- subset(dat, ok)					# 欠損値を持つケースを除く
	group <- group[ok]
	stopifnot(all(sapply(dat, is.factor), is.factor(group)))# データは全て factor であること
	nc <- nrow(dat)						# ケース数
	item <- ncol(dat)					# アイテム変数の数
	n <- as.vector(table(group))				# 各群のケース数
	ng <- length(n)						# 群の数
	vname <- colnames(dat)					# 変数名
	cname <- unlist(sapply(dat, levels))			# カテゴリー名
	gname <- levels(group)					# 群の名前
	dat <- data.frame(dat, group)
	dat[,1:(item+1)] <- lapply(dat, as.integer)
	cat <- sapply(dat, max)					# 各アイテム変数の最大値
	junjo <- c(0, cumsum(cat))[1:(item+1)]			# ダミー変数に変換するときに使う，アイテム変数の位置情報
	nobe2 <- sum(cat)					# 群変数も含めた延べカテゴリー数
	nobe <- nobe2-ng					# アイテム変数だけの延べカテゴリー数
	dat2 <- t(apply(dat, 1,					# 群変数も含めて，ダミー変数に変換する
			function(obs)
			{
				zeros <- numeric(nobe2)
				zeros[junjo+obs] <- 1
				zeros
			}
		))
	a2 <- array(apply(dat2, 1, function(i) outer(i, i)),	# クロス集計表を作る準備
			dim=c(nobe2, nobe2, nc))
	x <- apply(a2, 1:2, sum)				# 3 次元配列に対して合計をとる
	
	pcros <- x[1:nobe, 1:nobe]				# アイテム変数同士のクロス集計
	gcros <- x[(nobe+1):nobe2, 1:nobe]			# 群変数とアイテム変数間のクロス集計
	
	w <- outer(n, n)/nc
	diag(w) <- 0
	grgr <- diag(n-n^2/nc)-w
	
	w <- diag(pcros)
	grpat <- gcros-outer(n, w)/nc
	pat <- pcros-outer(w, w)/nc
	select <- (junjo+1)[1:item]				# 冗長な行の添え字
	
	pat <- pat[-select, -select]				# 冗長性を排除した行列
	grpat <- grpat[-1, -select, drop=FALSE]			# 冗長性を排除した行列
	r <- grgr[-1, -1, drop=FALSE]				# 冗長性を排除した行列
	
	ndim <- ng-1						# 得られる解の数
	m <- ncol(pat)						# 行列の大きさ
	
	pat <- solve(pat)					# 逆行列を求める
	
	c <- grpat%*%pat%*% matrix(t(grpat), nrow=m, ncol=ndim)
	
	w <- geneig(c, r)					# 一般化固有値問題を解く
	val <- w$values						# 固有値
	vec <- w$vectors					# 固有ベクトル
	
	w <- t(t(pat%*%t(grpat)%*%vec)/sqrt(val))
	a <- matrix(0, nrow=nobe, ncol=ndim)			# カテゴリースコア
	ie <- 0
	for (j in 1:item) {
		is <- ie+1
		ie <- is+cat[j]-2
		offset <- junjo[j]+2
		a[offset:(offset+ie-is),] <- w[is:ie,]
	}
	w <- diag(pcros)*a
	for (j in 1:item) {
		is <- junjo[j]+1
		ie <- is+cat[j]-1
		s <- colSums(w[is:ie,, drop=FALSE])/nc
		a[is:ie,] <- t(t(a[is:ie,])-s)
	}
	a <- a/sqrt(diag(t(a) %*% pcros %*% a)/nc)
	sample.score <- dat2[,1:nobe] %*% a			# サンプルスコア
	centroid <- n*rbind(0, vec)				# 各群の重心
	centroid <- t(t(rbind(0, vec))-colSums(centroid)/nc)
	centroid <- t(t(centroid)/sqrt(colSums(centroid^2*n)/nc/val))

	item1 <- item + 1					# 外的基準との偏相関係数
	partial.cor <- matrix(0, item, ndim)
	for (l in 1:ndim) {
		pat <- matrix(0, item1, item1)
		pat[1, 1] <- sum(n*centroid[, l]^2)
		temp <- rowSums(a[, l]*t(gcros*centroid[, l]))
		pat[1, 2:(item+1)] <- mapply(function(is, ie) sum(temp[is:ie]), junjo[-item1]+1, junjo[-item1]+cat[-item1])
		temp <- outer(a[, l], a[, l])*pcros
		for (i in 1:item) {
			is <- junjo[i]+1
			ie <- is+cat[i]-1
			for (j in 1:i) {
				pat[j+1, i+1] <- sum(temp[is:ie, (junjo[j]+1):(junjo[j]+cat[j])])
			}
		}
		d <- diag(pat)
		pat <- pat/sqrt(outer(d, d))
		pat <- pat+t(pat)
		diag(pat) <- 1
		pat <- solve(pat)
		partial.cor[, l] <- -pat[2:item1, 1]/sqrt(pat[1, 1]*diag(pat)[-1])
	}

	dim(val) <- c(1, ndim)
	rownames(a) <- paste(rep(vname, cat[1:item]), cname, sep=".")
	rownames(partial.cor) <- vname
	rownames(centroid) <- gname
	rownames(sample.score) <- paste("#", 1:nc, sep="")
	rownames(val) <- "Eta"
	colnames(val) <- colnames(a) <- colnames(centroid) <- colnames(sample.score) <- colnames(partial.cor) <- paste("解", 1:ndim, sep="")
	return(structure(list(ndim=ndim, group=group, ng=ng, category.score=a, partial.cor=partial.cor, 
		centroid=centroid, eta=val, sample.score=sample.score), class="qt2"))
}
# print メソッド
print.qt2 <- function(	obj,					# qt2 が返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	cat("\nカテゴリー・スコア\n\n")
	print(round(obj$category.score, digits=digits))
}
# summary メソッド
summary.qt2 <- function(obj,					# qt2 が返すオブジェクト
			digits=5)				# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.qt2 <- function(obj,					# qt2 が返すオブジェクト
	i=1,							# カテゴリースコアを描画するときの対象次元，三群以上の判別図のときの横軸に対する次元
	j=2,							# 三群以上の判別図のときの縦軸に対する次元
	xlab=colnames(obj$category.score)[i],			# 三群以上の判別図の横軸のラベル
	ylab=colnames(obj$category.score)[j],			# 三群以上の判別図の縦軸のラベル
			pch=1:obj$ng,				# 三群以上の scatterplot を描く記号
			col=1:obj$ng,				# 三群以上の scatterplot の記号の色
			xpos="topright", ypos=NULL,		# 三群以上の scatterplot の凡例の位置
	which=c("boxplot", "barplot", "category.score"),	# 描画するグラフの種類の指定（三群以上の判別のときはこの指示によらない）
	nclass=20,						# barplot の階級数
	...)							# plot, boxplot, barplot に引き渡す引数
{
	which <- match.arg(which)
	if (which == "category.score") {
		cscore <- obj$category.score
		cname <- rev(rownames(cscore))
		cscore <- rev(cscore[, i])
		names(cscore) <- NULL
		barplot(cscore, horiz=TRUE, xlab=xlab, ...)
		text(0, 1.2*(1:length(cname)-0.5), cname, pos=ifelse(cscore > 0, 2, 4))
	}
	else {
		group <- obj$group
		if (obj$ndim > 1) {
			group.levels <- levels(group)
			sample.score <- obj$sample.score
			plot(sample.score[, i], sample.score[, j], xlab=xlab, ylab=ylab, col=col[as.integer(group)], pch=pch[as.integer(group)], ...)
			legend(x=xpos, y=ypos, legend=group.levels, col=col, pch=pch)
		}
		else {
			if (which == "boxplot") {		# boxplot
				plot(obj$sample.score[, i] ~ group, xlab="群", ylab="サンプル・スコア", ...)
			}
			else if (which == "barplot") { 		# barplot
				tbl <- table(group, cut(obj$sample.score,
						breaks=pretty(obj$sample.score, n=nclass)))
				barplot(tbl, beside=TRUE, legend=TRUE, xlab="サンプル・スコア", ...)
			}
		}
	}
}
# 数量化 III 類
qt3 <-function(dat)							# カテゴリーデータ行列
{
	dat <- subset(dat, complete.cases(dat))				# 欠損値を持つケースを除く
	stopifnot(all(sapply(dat, is.factor)))
	vname <- colnames(dat)
	cname <- unlist(sapply(dat, levels))				# カテゴリー名
	item <- ncol(dat)
	dat[, 1:item] <- lapply(dat, as.integer)
	if (any(dat >= 3)) {						# データに 3 以上の値があれば，アイテムデータ
		cat <- sapply(dat, max)
		dat <- make.dummy(dat)					# アイテムデータをカテゴリーデータに変換
		colnames(dat) <- paste(rep(vname, cat), cname, sep=".")
	}
	else {
		dat <- as.matrix(dat)-1					# カテゴリーデータは 0/1 データに
		colnames(dat) <- vname					# 名前を付けておく
	}
	nc <- nrow(dat)							# ケース数
	nobe <- ncol(dat)						# 延べカテゴリー数
	ncc <- rowSums(dat)						# 行和
	ni <- colSums(dat)						# 列和
	if (any(ncc == 0)) stop("反応が0のケースがあります")
	if (any(ni == 0)) stop("反応が0のカテゴリーがあります")
	fnorm <- sum(ni)
	a2 <- array(apply(dat, 1, function(i) outer(i, i)), dim=c(nobe, nobe, nc))
	x <- 0
	junk <- sapply(1:nc, function(i) x <<- x+a2[,,i]/ncc[i])
	result <- eigen(x/sqrt(outer(ni, ni)))				# 固有値・固有ベクトルを求める
	ne <- length(result$values[result$values > 1e-5])		# 解の個数
	eval <- result$values[2:ne]					# 固有値
	evec <- result$vectors*sqrt(fnorm/ni)				# 固有ベクトル
	corr <- sqrt(eval)						# 相関係数
	cat.score <- evec[,2:ne]					# カテゴリー・スコア
	smp.score <- dat%*%cat.score/outer(ncc, sqrt(eval))		# サンプル・スコア
	rownames(cat.score) <- colnames(dat)
	names(eval) <- names(corr) <- colnames(cat.score) <- colnames(smp.score) <- paste("解", 1:(ne-1), sep="")
	rownames(smp.score) <- paste("#", 1:nc, sep="")
	return(structure(list(Eigen.value=eval, Correlation.coefficient=corr,
		Category.score=cat.score, Sample.score=smp.score), class="qt3"))
}
# print メソッド
print.qt3 <- function(	obj,
			ax=length(obj$Eigen.value),
			digits=5)
{
	res <- rbind(ev=obj$Eigen.value, eta=obj$Correlation.coefficient)
	rownames(res) <- c("固有値", "相関係数")
	print(res[, 1:ax], digits=digits)
	cat("\nカテゴリースコア\n\n")
	print(round(obj$Category.score[, 1:ax, drop=FALSE], digits=digits))
}
# summary メソッド
summary.qt3 <- function(obj,
			digits=5)
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.qt3 <- function(	obj,						# qt3 が返すオブジェクト
			axis.no=1:2,					# 横軸と縦軸にプロットする解の番号
			which=c("category.score", "sample.score"),	# カテゴリースコアを描くかサンプルスコアを描くかの指定
			xlab=NULL, ylab=NULL,				# 横軸，縦軸の名前
			axis=TRUE,					# 座標軸を描くかどうかの指定
			label.cex=0.7,					# 描画点につけるラベルの文字の大きさ
			...)						# plot 関数に引き渡す引数
{
	if (length(obj$Eigen.value) > 1) {
		which <- match.arg(which)
		if (which == "category.score") {
			dat <- obj$Category.score
		}
		else {
			dat <- obj$Sample.score
		}
		ax1 <- axis.no[1]
		ax2 <- axis.no[2]
		aname <- names(obj$Eigen.value)
		if (is.null(xlab)) {
			xlab <- aname[ax1]
		}
		if (is.null(ylab)) {
			ylab <- aname[ax2]
		}
		x <- dat[, ax1]
		y <- dat[, ax2]
		plot(x, y, xlab=xlab, ylab=ylab, ...)
		par(xpd=TRUE)
			labels <- rownames(dat)
		if (label.cex != 0) {
			text(x, y, labels, cex=label.cex, pos=1)
		}
		par(xpd=FALSE)
		if (axis) {
			abline(h=0, v=0)
		}
	}
	else {
		warning("解が一次元なので，二次元表示はできません。")
	}
}
# 数量化 IV 類
qt4 <- function(s)					# 類似度行列
{
	if (is.data.frame(s)) {
		s <- as.matrix(s)
	}
	n <- nrow(s)					# 行列のサイズ
	object.names <- colnames(s)
	if (is.null(object.names)) {
		object.names <- paste("対象", 1:n, sep="")
	}
	h <- s+t(s)					# 転置行列との和を計算して，
	diag(h) <- 0					# 対角要素を 0 にする
	diag(h) <- -rowSums(h)				# 行和を求めて新たな対角要素とする
	res <- eigen(h)					# 固有値固有ベクトルを求める
	values <- res$values[res$values > 1e-5]		# 固有値が 0.00001 以上のものを解とする
	ax <- length(values)				# 解の個数
	vectors <- res$vectors[,1:ax]			# 固有ベクトル
	colnames(vectors) <- names(values) <- paste("解", 1:ax, sep="") # 名前を付ける
	rownames(vectors) <- object.names		# 名前を付ける
	return(structure(list(ax=ax, n=n, values=values, vectors=vectors), class="qt4"))
}
# print メソッド
print.qt4 <- function(	res,				# princo が返すオブジェクト
			ax=res$ax,			# 何次元までの解を出力するか
			digits=5)			# 表示桁数
{
	ax <- min(ax, res$ax)
	val <- res$values
	val2 <- val/sum(val)
	val <- rbind(val, val2, cumsum(val2))
	rownames(val) <- c("固有値", "寄与率", "累積寄与率")
	print(round(val[, 1:ax], digits=digits))
	cat("\nベクトル\n\n")
	print(round(res$vectors[, 1:ax], digits=digits))
}
# plot メソッド
plot.qt4 <- function(	res,				# princo が返すオブジェクト
			text.cex=0.7,			# ラベルのフォントの大きさ
			...)				# plot への引数
{
	if (res$ax >= 2) {				# 二次元以上の解が得られたら，
		plot(res$vectors[,1:2], ...)		# 二次元の図を描く
		abline(v=0, h=0)
		old <- par(xpd=TRUE)
		text(res$vectors[,1], res$vectors[,2], rownames(res$vectors), pos=4, offset=.2, cex=text.cex)
		par(old)
	}
	else {
		warning("解が一次元なので，二次元配置図は描けません。")
	}
}
# 二次の判別関数
quad.disc <- function(	data,						# 説明変数データ行列
			group,						# グループを表すベクトル
			func.name=c("solve", "ginv"))			# 逆行列を求める関数
{
	inverse <- if (match.arg(func.name) == "solve") solve else { library(MASS); ginv}
	data <- as.data.frame(data)
	if (is.null(colnames(data))) {
		colnames(data) <- paste("Var", 1:p, sep="")
	}
	vname <- colnames(data)
	group <- as.factor(as.matrix(group))				# 群を表す変数は factor にする
	OK <- complete.cases(data, group)				# 欠損値を持つケースを除く
	data <- as.matrix(data[OK,])
	group <- group[OK]
	p <- ncol(data)							# 説明変数の個数
	n <- nrow(data)							# データの個数
	n.i <- table(group)						# 各群の例数
	g.name <- names(n.i)						# 群の名前
	k <- length(n.i)						# 群の個数
	group.means <- matrix(unlist(by(data, group, colMeans)), p)	# 各群・各変数の平均
	vars <- array(unlist(by(data, group, var)), c(p, p, k))		# 各群の分散・共分散行列
	inv.vars <- array(apply(vars, 3, inverse), c(p, p, k))		# 各群の分散・共分散行列の逆行列

	scores <- sapply(1:k, function(i) {				# 各ケースの各群からの距離
		temp <- t(data)-group.means[,i];
		sapply(1:n, function (j) temp[,j]%*%inv.vars[,,i]%*%temp[,j])
		}
	)
	p.values <- pchisq(scores, p, lower.tail=FALSE)			# 各ケースが各群に属する確率
	prediction <- as.factor(g.name[apply(p.values, 1, order)[k,]])	# どの群に属するか判別
	correct <- ifelse(prediction == group, TRUE, FALSE)		# 判別が正しいか？
	correct.table <- table(group, prediction)			# 判別結果概括表
	correct.rate <- sum(diag(correct.table))/n*100			# 正判別率
# add names
	colnames(group.means) <- colnames(scores) <- colnames(p.values) <- dimnames(vars)[[3]] <- dimnames(inv.vars)[[3]] <- g.name
	colnames(vars) <- rownames(vars) <- colnames(inv.vars) <- rownames(inv.vars) <- rownames(group.means) <- vname
	rownames(scores) <- rownames(p.values) <- paste("case", 1:n)
	return(structure(list(group.means=group.means, vars=vars, inv.vars=inv.vars, scores=scores, p.values=p.values, prediction=prediction, correct=correct, correct.table=correct.table, correct.rate=correct.rate, group=group, ngroup=k), class="quad.disc"))
}
# print メソッド
print.quad.disc <- function(	obj,					# quad.disc 関数が返すオブジェクト
				digits=5)				# 結果の表示桁数
{
	cat("\n判別結果\n\n")
	print(obj$correct.table)
	cat(sprintf("\n正判別率 = %.1f %%\n", obj$correct.rate))
}
# summary メソッド							# すべての結果を表示する
summary.quad.disc <- function(	obj,					# quad.disc が返すオブジェクト
				digits=5)				# 結果の表示桁数
{
	print.default(obj, digits=digits)
}
# plot メソッド
plot.quad.disc <- function(	obj, 					# quad.disc 関数が返すオブジェクトの
				which=c("boxplot", "barplot", "scatterplot"),	# 箱髭図か棒グラフか散布図かの選択
				nclass=20,				# barplot の場合のおよその階級数
				pch=1:obj$ngroup,			# scatterplot を描く記号
				col=1:obj$ngroup,			# scatterplot の記号の色
				xpos="topright", ypos=NULL,		# scatterplot の凡例の位置
				...)					# boxplot, barplot, plot に引き渡す引数
{
	if (obj$ngroup == 2) {
		group <- obj$group
		score <- obj$score[,1]-obj$score[,2]
		which <- match.arg(which)
		if (which == "boxplot") {				# boxplot
			plot(score ~ group, xlab="群", ylab="二乗距離の差", ...)
		}
		else if (which == "barplot") { 				# barplot
			tbl <- table(group, cut(score,
					breaks=pretty(score, n=nclass)))
			barplot(tbl, beside=TRUE, legend=TRUE, xlab="二乗距離の差", ...)
		}
		else {							# scatterplot 各群の重心までの二乗距離
			group <- obj$group
			group.levels <- levels(group)
			score1 <- obj$score[,1]
			score2 <- obj$score[,2]
			max1 <- max(score1)
			max2 <- max(score2)
			max0 <- max(max1, max2)
			plot(score1, score2, col=col[as.integer(group)], pch=col[as.integer(group)],
				xlim=c(0, max0), xlab=paste(group.levels[1], "の重心への二乗距離"),
				ylim=c(0, max0), ylab=paste(group.levels[2], "の重心への二乗距離"), asp=1, ...)
			abline(0, 1, lty=2)
			text(max1, max2/2, paste(group.levels[2], "に判別される領域"), pos=2)
			text(0, max2+strheight("H")*1.5, paste(group.levels[1], "に判別される領域"), pos=4)
			legend(x=xpos, y=ypos, legend=group.levels, col=col, pch=pch)
		}
	}
	else {
		warning("3群以上の場合にはグラフ表示は用意されていません")
	}
}
radar <- function(	df,		# データフレームまたはデータ行列
			max=NULL,	# 変数値の上限（変数単位ならベクトルで，全部一緒なら定数で）NULL なら計算
			min=NULL,	# 変数値の下限（変数単位ならベクトルで，全部一緒なら定数で）NULL なら計算
			z.score=TRUE,	# データを正規化してプロット
			col=2,		# 線の色（データ単位ならベクトルで，全部一緒なら定数で）
			lty=1,		# 線の種類（データ単位ならベクトルで，全部一緒なら定数で）
			title="")	# グラフのタイトル
{
	draw.net <- function(x, border, lty)			# 1個のデータの描画
	{
		scale <- (x-min)/(max-min)			# 描画に使う値に変換
		if (0 < scale && scale <= 1.2) {		# 範囲外のときのみ描画（z.score=TRUE のときのみ有効）
			polygon(sine*scale, cosine*scale,	# 多角形を描画
				border=border, lty=lty)
		}	
	}

	if (!is.data.frame(df) && !is.matrix(df)) {		# オブジェクトの種類
		stop("行列かデータフレームであるべし！")
	}
	m <- ncol(df)						# 変数の個数
	if (m < 3) {
		stop("3変数以上であるべし！")
	}
	plot(c(-1.2,1.2),c(-1.2,1.2), type="n", axes=FALSE,	# 描画の枠組み
	        xlab="", ylab="", main=title, asp=1)
	theta <- pi*(0.5-0:(m-1)*2/m)				# 90 度から右回りで測る角度
	sine <- cos(theta)					# 正弦
	cosine <- sin(theta)					# 余弦
	if (z.score) {						# 正規化してプロットするとき
		df <- scale(df)					# 変数ごとに正規化
		max <- max(df)					# 変数ごとの最大値
		min <- min(df)					# 変数ごとの最小値
		w <- (max-min)*0.1				# マージン
		max <- max+w					# マージンを拡大
		min <- min-w					# マージンを拡大
		sapply(-3:3, function(i)			# 1ごとに目盛り枠を描画
			draw.net(i, border="gray", lty=ifelse(i==0, 1, 3)))
	}
	else {							# 生データのままプロット
		if (is.null(max)) {				# 最大値の指定がないとき
			max <- apply(df, 2, max)		# 最大値を求める
		}
 		if (is.null(min)) {				# 最小値の指定がないとき
 			min <- apply(df, 2, min)		# 最小値を求める
 		}
		w <- (max-min)*0.1				# マージン
		max <- max+w					# マージンを拡大
		min <- min-w					# マージンを拡大
		sapply(1:5, function(i)				# 目盛り枠を描画
			polygon(sine*i/5, cosine*i/5, lty=3, border="gray"))
	}
	arrows(0, 0, sine*1.1, cosine*1.1, length=0,		# 軸を描画
		col="gray")
	half <- floor((m+1)/2)					# 右半分と左半分を分ける
	text(sine*1.2, cosine*1.2, labels=colnames(df),		# ちょっと手を掛ける
		pos=rep(c(4, 2), c(half, m-half)))
	n <- nrow(df)						# サンプルサイズ
	if (length(col) == 1) {					# 要素が 1 個のとき
		col <- rep(col, n)				# ベクトルに拡大
	}
	if (length(lty) == 1) {					# 要素が1個のとき
		lty <- rep(lty, n)				# ベクトルに拡大
	}
	junk <- sapply(1:n, function(i)				# 描画
		draw.net(df[i,], col[i], lty[i]))
}
# 乱塊法
randblk <- function(dat)				# データ行列
{
	dat <- as.matrix(dat)
	dat <- subset(dat, complete.cases(dat))		# 欠損値を持つケースを除く
	nr <- nrow(dat)					# ケース数
	nc <- ncol(dat)					# 変数の個数（処理・条件の個数）
	gmean <- mean(dat)				# 全平均値
	tmean <- colMeans(dat)				# 処理・条件ごとの平均値
	rmean <- rowMeans(dat)				# 繰り返しの平均値
	
	SSt <- sum((dat-gmean)^2)
	SSc <- nr*sum((tmean-gmean)^2)
	SSr <- nc*sum((rmean-gmean)^2)
	SSe <- SSt-SSr-SSc
	SS <- c(SSc, SSr, SSe, SSt)

	dfc <- nc-1
	dfr <- nr-1
	dfe <- dfc*dfr
	dft <- nc*nr-1
	df <- c(dfc, dfr, dfe, dft)

	MS <- SS/df

	Fs <- MS/MS[3]
	Ps <- pf(Fs, df, dfe, lower.tail=FALSE)
	Fs[3:4] <- Ps[3:4] <- NA

	result <- cbind(SS, df, MS, Fs, Ps)
	rownames(result) <- c("Treatment", "Replication", "Residual", "Total")
	colnames(result) <- c("SS", "d.f.", "MS", "F value", "P value")
	return(structure(list(result=result), class="randblk"))
}
# print メソッド					分散分析表を表示
print.randblk <- function(	obj,			# randblk 関数が返すオブジェクト
				digits=5)		# 表示桁数
{
	print(obj$result, na.print="", digits=digits)
}
# 連続変数データをカテゴリーデータに変換
Recode <- function(	x,							# 連続変数データベクトル
			arg1, arg2)						# 二通りの意味を持つ
{
	if (length(arg1) == 1 && length(arg2) == 1) {				# arg1: 最小値が含まれる区間の左限界値
		return(factor(floor((x-arg1)/arg2)*arg2+arg1))			# arg2: 区間幅
	}
	else if (length(arg1) == length(arg2)-1) {				# arg1: カッティングポイント
		cut(x, breaks=c(-Inf, arg1, Inf), right=FALSE, labels=arg2)	# arg2: 区切られた区間のラベル
	}
	else {
		stop("区切り値の個数は，ラベルの数よりちょうど 1 だけ小さいはずです")
	}
}
# カテゴリーデータの再カテゴリー化
recode2 <- function(	x,				# 連続変数データベクトル
			y)				# 再カテゴリー化を定義する n×3 行列
							# i 行の 1 列目 2 列目を y[i,1], y[i,2] としたとき，
							# y[i,1] ≦ x ≦ y[i,2] である x を y[i,3] に置き換える
{
	if (ncol(y) != 3) {
		stop("y は3列でなければいけない")
	}
	n <- nrow(y)
	for (i in 1:n) {
		for (j in 1:n) {
			if (i != j && ((y[i,1] <= y[j,1] && y[j,1] <= y[i,2]) || (y[i,1] <= y[j,2] && y[j,2] <= y[i,2]))) {
				stop("区間情報に重複があります")
			}
		}
	}
	sapply(x, function(z) ifelse(is.na(z), NA, y[,3][y[,1] <= z & z <= y[,2]]))
}
# 相対危険度（対応のない場合）とその信頼限界を求める
relative.risk <- function(	a,			# 対象群・所見あり
				b,			# 対象群・所見なし
				c,			# 対照群・所見あり
				d)			# 対照群・所見なし
{
	cl <- function(x)
	{
		exp(log(rr)+c(1, -1)*qnorm(x)*sqrt(b/a/(a+b)+d/c/(c+d)))
	}
	rr <- a*(c+d)/c/(a+b)				# 相対危険度
	conf <- rbind(cl90=cl(0.05), cl95=cl(0.025), cl99=cl(0.005), cl999=cl(0.0005))
	colnames(conf) <- paste(c("下側","上側"), "信頼限界値", sep="")
	rownames(conf) <- paste(c(90, 95, 99, 99.9), "%信頼区間", sep="")
	list(rr=rr, conf=conf)
}
# 相対危険度（対応のある場合）とその信頼限界を求める
relative.risk2 <- function(	b,			# 対照群要因あり，症例群要因なし
				c)			# 対照群要因なし，症例群要因あり
{
	cl <- function(x)
	{
		exp(log(rr)+c(1, -1)*qnorm(x)*sqrt(1/b+1/c))
	}
	rr <- c/b					# 相対危険度
	conf <- rbind(cl90=cl(0.05), cl95=cl(0.025), cl99=cl(0.005), cl999=cl(0.0005))
	colnames(conf) <- paste(c("下側","上側"), "信頼限界値", sep="")
	rownames(conf) <- paste(c(90, 95, 99, 99.9), "%信頼区間", sep="")
	list(rr=rr, conf=conf)
}
# 抵抗直線を描く
resistant.line <- function(	x,					# 独立変数（横軸に取る）
				y,					# 従属変数（縦軸に取る）
				no.iteration=FALSE,			# 収束計算をするかどうか
				n=1)					# ブートストラップ法で信頼区間を求めるときの回数
{
	resistant.line0 <- function(x, y)
	{
		OK <- complete.cases(x, y)				# 欠損値を持つケースを除く
		x <- x[OK]
		y <- y[OK]

		n <- length(x)						# データの組数
		o <- order(x)						# x の順序づけ
		x <- x[o]						# x を昇順に並べ替える
		y <- y[o]						# x，y を対応づけて並べ替える
		m <- floor(n/3)						# データを 3 分割するときの個数
		n1 <- 1:m						# 最初の 1/3 に対する添え字
		n2 <- (m+1): (n-m)					# 次の 1/3 に対する添え字
		n3 <- (n-m+1):n						# 最後の 1/3 に対する添え字
	
		m.x1 <- median(x[n1])					# 最初の 1/3 の中央値
		m.x2 <- median(x[n2])					# 次の 1/3 の中央値
		m.x3 <- median(x[n3])					# 最後の 1/3 の中央値
	
		a <- b <- 0						# 傾きと切片の初期値
		repeat {						# 収束計算
			m.y1 <- median(y[n1])				# 最初の 1/3 の中央値
			m.y2 <- median(y[n2])				# 次の 1/3 の中央値
			m.y3 <- median(y[n3])				# 最後の 1/3 の中央値
			d.a <- (m.y3-m.y1)/(m.x3-m.x1)			# 傾きの修正量
			d.b <- (m.y1+m.y2+m.y3-d.a*(m.x1+m.x2+m.x3))/3	# 切片の修正量
			a <- a+d.a					# 傾きの修正
			b <- b+d.b					# 切片の修正
			if (no.iteration || abs(d.a/a) < 1e-7) break	# 収束計算をしない場合，収束した場合はループ脱出
			y <- y-(d.a*x+d.b)				# 抵抗直線で説明できる部分を取り除く
		}
		return(c(b, a))						# 結果を返す
	}
	Driver <- function(x, y)					# ブートストラップ法のためのドライバー
	{
		n <- length(x)
		suffix <- sample(n, n, replace=TRUE)			# リサンプリング
		return(resistant.line0(x[suffix], y[suffix]))		# リサンプリングしたデータについてパラメータを推定
	}
	names.xy <- c(deparse(substitute(x)), deparse(substitute(y)))	# 変数名を控えておく
	ans <- list(coefficients=resistant.line0(x, y),			# 引数に対してパラメータを推定する
		    names.xy=names.xy, x=x, y=y)
	if (n > 1) {
		ans2 <- replicate(n, Driver(x, y))			# ブートストラップを n 回実行
		ans <- append(ans, list(intercepts=ans2[1,], slopes=ans2[2,]))
	}
	class(ans) <- "resistant.line"					# print, plot メソッドのためにクラス名をつけておく
	return(ans)
}
# print メソッド
print.resistant.line <- function(	obj,				# "resistant.line" オブジェクト
					digits=5,			# 表示桁数
					sig=0.95)			# 信頼度
{
	if (length(obj) == 4) {
		cat("Intercept:", round(obj$coefficients[1], digits),
		    "    Slope:", round(obj$coefficients[2], digits), "\n")
	}
	else {
		alpha <- (1-sig)/2
		LCL <- c(quantile(obj$intercepts, alpha), quantile(obj$slopes, alpha))
		UCL <- c(quantile(obj$intercepts, 1-alpha), quantile(obj$slopes, 1-alpha))
		ans <- data.frame(obj$coefficients, LCL, UCL)
		dimnames(ans) <- list(c("Intercept:", "Slope:"),
				      c("Estimate", paste(c(alpha, 1-alpha), "%", sep="")))
		print(ans, digits=digits)		
	}
}
# plot メソッド
plot.resistant.line <- function(obj,					# "resistant.line" オブジェクト
			posx="topleft", posy=NULL,			# legend 関数のための位置引数
			xlab=obj$names.xy[1], ylab=obj$names.xy[2],	# 軸の名前
			hist=FALSE,					# ヒストグラムを描くとき TRUE にする
			...)						# その他の任意の plot 関数の引数
{
	if (hist && length(obj) == 6) {					# ブートストラップの結果を，hist=TRUE のときに，ヒストグラムで表示する
		layout(matrix(1:2, 2))
		hist(obj$intercepts, xlab="Intercept", main="", right=FALSE)
		hist(obj$slopes, xlab="Slope", main="", right=FALSE)
		layout(1)
	}
	else {								# 散布図と Deming 法の回帰直線と直線回帰式を表示する
		plot(obj$x, obj$y, xlab=xlab, ylab=ylab, ...)
		abline(obj$coefficients)
		abline(lm(obj$y~obj$x), lty=2, col=2)
		legend(posx, posy, legend=c("resistant line", "linear regression"), lty=1:2, col=1:2)
	}
}
# 順位データを双対尺度法で分析する
ro.dual <- function(F)						# 順位データ
{
	F <- data.matrix(F)						# データフレームも行列にする
	N <- nrow(F)							# 評価者の数
	if (is.null(rownames(F))) {					# 行名（評価者名）がないとき，
		row.names <-  paste("Row", 1:N, sep="-")		# 行名の補完
	}
	n <- ncol(F)							# 評価対象の数
	if (is.null(colnames(F))) {					# 列名（評価対象名）がないとき，
		col.names <- paste("Col", 1:n, sep="-")		# 列名の補完
	}
	E <- n+1-2*F
	Hn <- t(E)%*%E/(N*n*(n-1)^2)
	ans <- eigen(Hn)						# 固有値・固有ベクトルを求める
	ne <- nrow(Hn)-1						# 有効な固有値・固有ベクトルの個数
	eta2 <- ans$values[1:ne]					# 固有値（相関比の二乗）
	eta <- sqrt(eta2)						# 相関比
	contribution <- eta2/sum(ans$values[1:ne])*100		# 寄与率
	cumcont <- cumsum(contribution)				# 累積寄与率
	result <- rbind(eta2, eta, contribution, cumcont)	# 結果
	dimnames(result) <- list(c("eta square", "correlation", "contribution", "cumulative contribution"),
				 paste("Axis", 1:ne, sep="-"))
	W <- ans$vectors[, 1:ne, drop=FALSE]			# 固有ベクトル
	col.score <- W*sqrt(n)					# 列スコア
	col.score2 <- t(t(col.score)*eta)				# 相関比で重み付けした列スコア
	row.score2 <- t(t(E%*%W/sqrt(n)/(n-1)))			# 相関比で重み付けした行スコア
	row.score <- t(t(row.score2)/eta)				# 行スコア
	colnames(col.score) <- colnames(row.score) <- colnames(result)
	rownames(col.score) <- col.names
	rownames(row.score) <- row.names
	dimnames(col.score2) <- dimnames(col.score)
	dimnames(row.score2) <- dimnames(row.score)
	result <- list(	result=result,
			row.score=row.score,
			col.score=col.score, 
			row.score.weighted=row.score2, 
			col.score.weighted=col.score2)
	class(result) <- "dual"					# summary, plot メソッドがある
	return(result)
}	
# 母相関係数が 0 以外の特定の値であるかどうかの検定
rtest <- function(	n,			# 標本の大きさ
			r,			# 標本相関係数
			rho)			# 母相関係数
{
	method <- "母相関係数が 0 以外の特定の値であるかどうかの検定"
	data.name <- sprintf("n = %s, r = %s, rho = %s", n, r, rho)
	z <- abs(atanh(r)-atanh(rho))*sqrt(n-3)	# 検定統計量
	p <- pnorm(z, lower.tail=FALSE)*2	# P 値
	return(structure(list(statistic=c(z=z),	# 結果をまとめて返す
		p.value=p, method=method, data.name=data.name),
		class="htest"))
}
# 二群の平均値の差（両側検定）を行うときに必要な各群あたりのサンプルサイズを求める
sample.size <- function(alpha,					# 有意水準
			powd,					# 検出力
			esize)					# 効果量
{

	gcf <- function(a, x)
	{
		ITMAX <- 100
		EPS <- 3e-7
		FPMIN <- 1e-30
	
		b <- x+1-a
		c <- 1/FPMIN
		d <- 1/b
		h <- d
		for (i in 1:ITMAX) {
			an <- -i*(i-a)
			b <- b+ 2
			d <- an*d+b
			if (abs(d) < FPMIN) d <- FPMIN
			c <- b+an/c
			if (abs(c) < FPMIN) c <- FPMIN
			d <- 1/d
			del <- d*c
			h <- h*del
			if (abs(del-1) < EPS) {
				return(exp(-x+a*log(x)-lgamma(a))*h)
			}
		}
		stop("error")
	}

	gser <- function(a, x)					# 不完全ガンマ関数
	{
		ITMAX <- 100
		EPS <- 3e-7
	
		if (x == 0) {
			return(0)
		}
		else if (x > 0) {
			ap <- a
			del <- sum <- 1/a
			for (n in 1:ITMAX) {
				ap <- ap+1
				del <- del*x/ap
				sum <- sum+del
				if (abs(del) < abs(sum)*EPS) {
					return(sum*exp(-x+a*log(x)-lgamma(a)))
				}
			}
		}
		stop("error")
	}

	gammp <- function(a, x)					# pgamma(x, a)
	{
		ifelse(x < a+1, gser(a, x), 1-gcf(a, x))
	}

	erff <- function(x)					# 1-2*pnorm(x*sqrt(2), lower.tail=FALSE)
	{
		ifelse(x < 0, -gammp(0.5, x*x), gammp(0.5, x*x))
	}

	betacf <- function(a, b, x)
	{
		ITMAX <- 100
		EPS <- 3e-7
		FPMIN <- 1e-30
	
		qab <- a+b
		qap <- a+1
		qam <- a-1
		c <- 1
		d <- 1-qab*x/qap
		if (abs(d) < FPMIN) d <- FPMIN
		d <- 1/d
		h <- d
		for (m in 1:ITMAX) {
			m2 <- 2*m
			aa <- m*(b-m)*x/((qam+m2)*(a+m2))
			d <- 1+aa*d
			if (abs(d) < FPMIN) d <- FPMIN
			c <- 1+aa/c
			if (abs(c) < FPMIN) c <- FPMIN
			d <- 1/d
			h <- h*d*c
			aa <- -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
			d <- 1+aa*d
			if (abs(d) < FPMIN) d <- FPMIN
			c <- 1+aa/c
			if (abs(c) < FPMIN) c <- FPMIN
			d <- 1/d
			del <- d*c
			h <- h*del
			if (abs(del-1) < EPS) return(h)
		}
		stop("error")
	}

	betai <- function(a, b, x)				# 不完全ベータ関数
	{
		if (x < 0 || x > 1) stop("error")
		bt <- ifelse(x == 0 || x == 1, 0, exp(lbeta(a, b)+a*log(x)+b*log(1-x)))
		ifelse(x < (a+1)/(a+b+2), bt*betacf(a, b, x)/a, 1-bt*betacf(b, a, 1-x)/b)
	}

	sub1 <- function(n, esize, alpha)
	{
		df <- n-2
		t <- qt(alpha/2, df, lower.tail=FALSE)
		dd <- 1-0.25/df+1/(32*df*df)
		0.5*(1+erff((esize*sqrt(n)-sqrt(2)*t*dd) / (2*sqrt(1+t*t*(1-dd*dd)))))
	}

# 関数本体

	n <- 0
	powa <- 0
	INTV <- 200
	EPS <- 0.001
	dir <- -1

	while (powa <= powd) {
		n <- n+100
		powa <- sub1(n, esize, alpha)
	}

	while (abs((powa-powd)/powd) >= EPS) {
		INTV <- INTV*0.5
		n <- n+dir*INTV*0.5
		powa <- sub1(n, esize, alpha)
		dir <- ifelse(powa < powd, 1, -1)
	}
	return(n)
}
# 散布図を描き，棄却楕円（確率楕円），回帰直線，回帰直線の信頼限界帯，MA，RMA による回帰直線を描く
scatter <- function(	x,					# 独立変数（横軸にとる）
			y,					# 従属変数（縦軸にとる）
			ellipse=FALSE,				# 確率楕円を加えるときに TRUE *1
			lrl=FALSE,				# 回帰直線を加えるときに TRUE *1
			cb=FALSE,				# 回帰直線の信頼限界を加えるときに TRUE
			ma=FALSE,				# MA による回帰直線を加えるときに TRUE *1
			rma=FALSE,				# RMA による回帰直線を加えるときに TRUE *1
			alpha=0.05,				# 1-α が信頼度（信頼率）
			acc=2000,				# 確率楕円のなめらかさ
			xlab=NULL,				# x 軸名称
			ylab=NULL)				# y 軸名称
								# *1 TRUE/FALSE 以外に，長さ 3 の数値ベクトルで，
								# (線種, 太さ, 色) を指定できる。数値は lty, lwd, col を参照のこと
{								# 棄却楕円のときは，色は，border の色
	comp <- function(a)
	{
		if (length(a) == 1) a[2] <- a[3] <- 1
		if (length(a) == 2) a[3] <- 1
		return(a)
	}	
	
	ellipse.draw <- function(x, y, alpha, acc, ellipse, 	# 棄却楕円を描く
				 xlab, ylab)
	{
		ellipse <- comp(ellipse)
		vx <- var(x)
		vy <- var(y)
		vxy <- var(x, y)
		lambda <- eigen(var(cbind(x, y)))$values
		a <- sqrt(vxy^2/((lambda[2]-vx)^2+vxy^2))
		b <- (lambda[2]-vx)*a/vxy
		theta <- atan(a/b)
		k <- sqrt(-2*log(alpha))
		l1 <- sqrt(lambda[1])*k
		l2 <- sqrt(lambda[2])*k
		x2 <- seq(-l1, l1, length.out=acc)
		tmp <- 1-x2^2/l1^2
		y2 <- l2*sqrt(ifelse(tmp < 0, 0, tmp))
		x2 <- c(x2, rev(x2))
		y2 <- c(y2, -rev(y2))
		s0 <- sin(theta)
		c0 <- cos(theta)
		xx <- c0*x2+s0*y2+mean(x)
		yy <- -s0*x2+c0*y2+mean(y)
		rngx <- range(c(x, xx))
		rngy <- range(c(y, yy))
		plot(rngx, rngy, type="n", xlab=xlab, ylab=ylab)
		polygon(xx, yy, lty=ellipse[1], lwd=ellipse[2], border=ellipse[3])
	}

	conf.limit <- function(x, y, cb, alpha, lrl)		# 回帰直線と信頼限界帯を描く
	{
		lrl <- comp(lrl)
		n <- length(x)
		b <- var(x, y)/var(x)
		a <- mean(y)-b*mean(x)
		abline(a, b, lty=lrl[1], lwd=lrl[2], col=lrl[3])
		if (cb) {
			sx2 <- var(x)*(n-1)
			R <- par()$usr				# 横軸の範囲
			x1 <- seq(R[1], R[2], length.out=2000)
			y1 <- a+b*x1
			ta <- -qt(alpha/2, n-2)
			Ve <- (var(y)-var(x, y)^2/var(x))*(n-1)/(n-2)
			temp <- ta*sqrt(Ve)*sqrt(1/n+(x1-mean(x))^2/sx2)
			y2 <- y1-temp
			lines(x1, y2, lty="dotted")
			y2 <- y1+temp
			lines(x1, y2, lty="dotted")
			temp <- ta*sqrt(Ve)*sqrt(1+1/n+(x1-mean(x))^2/sx2)
			y2 <- y1-temp
			lines(x1, y2, lty="dashed")
			y2 <- y1+temp
			lines(x1, y2, lty="dashed")
		}
		cat("LS(least squares)--------------------\n")
		list(intercept=a, slope=b)
	}
	
	MA <- function(x, y, ma)				# Major Axis regression
	{
		ma <- comp(ma)
		s2 <- cov(cbind(x, y))
		b <- s2[1, 2]/(eigen(s2)$values[1]-s2[2, 2])
		a <- mean(y)-b*mean(x)
		abline(a, b, lty=ma[1], lwd=ma[2], col=ma[3])
		cat("MA(major axis)--------------------\n")
		list(intercept=a, slope=b)
	}
	
	RMA <- function(x, y, rma)					# Reduced Major Axis regression
	{
		rma <- comp(rma)
		b <- sign(cor(x, y))*sqrt(var(y)/var(x))
		a <- mean(y)-b*mean(x)
		abline(a, b, lty=rma[1], lwd=rma[2], col=rma[3])
		cat("RMA(reduced major axis)--------------------\n")
		list(intercept=a, slope=b)
	}

	if (is.null(xlab)) xlab <- deparse(substitute(x))	# x 軸名称
	if (is.null(ylab)) ylab <- deparse(substitute(y))	# y 軸名称
	OK <- complete.cases(x, y)				# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	if (ellipse[1]) {					# データポイントをマークして棄却楕円を描く
		ellipse.draw(x, y, alpha, acc, ellipse, xlab, ylab)
		points(x, y)
	}
	else {							# 散布図のみ
		plot(x, y, xlab=xlab, ylab=ylab)
	}
	if (lrl[1]) {						# 回帰直線と信頼限界帯
		print(conf.limit(x, y, cb, alpha, lrl))
	}
	if (ma[1]) {						# Major Axis regression
		print(MA(x, y, ma))
	}
	if (rma[1]) {						# Reduced Major Axis regression
		print(RMA(x, y, rma))
	}
}
# シェッフェの方法による線形比較
scheffe <- function(	n,			# 各群のデータ個数のベクトル
			m,			# 各群の平均値のベクトル
			u,			# 各群の不偏分散のベクトル
			g1,			# 第一グループの指定
			g2,			# 第二グループの指定
			conf.level=0.95)	# 線形比較の信頼区間を計算する信頼率
{
	stopifnot(length(n) == length(m), length(m) == length(u), n > 1, u > 0, floor(n) == n, floor(g1) == g1, floor(g2) == g2)
	method <- "シェッフェの方法による線形比較"
	data.name <- paste(deparse(substitute(g1)), "and", deparse(substitute(g2)))
	ng <- length(n)				# 群の数
	k1 <- ng-1
	nc <- sum(n)				# 全体のデータ数
	dfw <- nc-ng				# 群内平方和の自由度
	Vw <- sum(u*(n-1))/dfw			# 群内平均平方（群内不偏分散）
	n1 <- length(g1)			# 第一グループにまとめる群数
	n2 <- length(g2)			# 第二グループにまとめる群数
	g0 <- (1:ng)[-c(g1, g2)]		# どちらのグループにも含まれない群の番号
	n0 <- ng-n1-n2				# どちらのグループにも含まれない群の数
	weight <- rep(c(1/n1, -1/n2, 0),	# 重み（合計すると 0 になる）
			c(n1, n2, n0))[order(c(g1, g2, g0))]
	theta <- sum(weight*m)			# θ推定量
	Vtheta <- Vw*sum(weight^2/n)		# θの分散
	conf.int <- theta-c(1, -1)*		# θの信頼区間
		sqrt(k1*qf(1-conf.level, k1, dfw, lower.tail=FALSE)*Vtheta)
	attr(conf.int, "conf.level") <- conf.level
	F0 <- theta^2/k1/Vtheta			# F 値
	p <- pf(F0, k1, dfw, lower.tail=FALSE)	# P 値
	return(structure(list(statistic=c(theta=theta, "V(theta)"=Vtheta, F=F0),
		parameter=c(df1=k1, df2=dfw), p.value=p, conf.int=conf.int,
		method=method, data.name=data.name, contrast=list(g1, g2)), class="htest"))
}
sdis <- function(	data,					# 説明変数だけのデータフレーム
			group,					# 群を表す変数（ベクトルではなく，1 列のデータフレームとして引用するほうがよい）
			stepwise=TRUE,				# ステップワイズ変数選択をする
			P.in=0.05,				# Pin
			P.out=0.05,				# Pout （Pout ≧ Pin のこと）
			predict=FALSE,				# 各ケースの予測値を出力する
			verbose=TRUE)				# ステップワイズ変数選択の途中結果を出力する
{
	step.out <- function(isw)				# 変数選択の途中結果を出力
	{
		step <<- step+1					# ステップ数
		ncase.k <- ncase-ng
		if (isw != 0 && verbose) {
			cat(sprintf("\n***** ステップ %i *****   ", step))
			cat(sprintf("%s変数: %s\n", c("編入", "除去")[isw], vname[ip]))
		}
		lxi <- lx[1:ni]
		lxi2 <- cbind(lxi, lxi)
		a <- matrix(0, ni, ng)
		a0 <- numeric(ng)
		for (g in 1:ng) {
			a[, g] <- -(w[lxi, lxi]%*%Mean[lxi, g])*2*ncase.k
			a0[g] <- Mean[lxi, g]%*%w[lxi, lxi]%*%Mean[lxi, g]*ncase.k
		}
		idf1 <- ng-1
		idf2 <- ncase-(ng-1)-ni
		temp <- idf2/idf1
		f <- t[lxi2]/w[lxi2]				# 偏 F 値
		f <- temp*(1-f)/f
		P <- pf(f, idf1, idf2, lower.tail=FALSE)	# P 値
		rownames(a) <- c(vname[lxi])
		result2 <- data.frame(rbind(a, a0), f=c(f, NA), p=c(format.pval(P, 3, 1e-3), NA))
		dimnames(result2) <- list(c(vname[lxi], "定数項"), c(grp.lab, "偏F値", "P値"))
		class(result2) <- c("sdis", "data.frame")	# print.sdis を使うための設定
		if (verbose) {
			cat("\n***** 分類関数 *****\n")
			print(result2)
		}
		alp <- ng-1
		b <- ncase-1-0.5*(ni+ng)
		qa <- ni^2+alp^2
		c <- 1
		if (qa != 5) {
			c <- sqrt((ni^2*alp^2-4)/(qa-5))
		}
		df1 <- ni*alp					# 第1自由度
		df2 <- b*c+1-0.5*ni*alp				# 第2自由度
		wl <- detw/dett					# ウィルクスの Λ
		cl <- exp(log(wl)/c)
		f <- df2*(1-cl)/(df1*cl)			# 等価な F 値
		p <- pf(f, df1, df2, lower.tail=FALSE)		# P 値
		if (verbose) {
			cat(sprintf("ウィルクスのΛ: %.5g\n", wl))
			cat(sprintf("等価なＦ値:　　 %.5g\n", f))
			cat(sprintf("自由度:　　　　 (%i, %.2f)\n", df1, df2))
			cat(sprintf("Ｐ値:　　　　　 %s\n", format.pval(p, 3, 1e-3)))
		}
		return(result2)
	}
	
	fmax <- function()					# モデルに取り入れる変数の探索
	{
		kouho <- 1:p
		if (ni > 0) {
			kouho <- (1:p)[-lx[1:ni]]
		}
		kouho <- cbind(kouho, kouho)
		temp <- w[kouho]/t[kouho]
		temp <- (1-temp)/temp
		ip <- which.max(temp)
		return(c(temp[ip], kouho[ip]))
	}
	
	fmin <- function()					# モデルから捨てる変数の探索
	{
		kouho <- cbind(lx[1:ni], lx[1:ni])
		temp <- t[kouho]/w[kouho]
		temp <- (1-temp)/temp
		ip <- which.min(temp)
		return(c(temp[ip], lx[ip]))
	}
	
	sweep.sdis <- function(r, det)				# 掃き出し法
	{
		ap <- r[ip, ip]
		if (abs(ap) <= EPSINV) {
			stop("正規方程式の係数行列が特異行列です")
		}
		det <- det*ap
		for (i in 1:p) {
			if (i != ip) {
				temp <- r[ip, i]/ap
				for (j in 1:p) {
					if (j != ip) {
						r[j, i] <- r[j, i]-r[j, ip]*temp
					}
				}
			}
		}
		r[,ip] <- r[,ip]/ap
		r[ip,] <- -r[ip,]/ap
		r[ip, ip] <- 1/ap
		return(list(r=r, det=det))
	}
	
	discriminant.function <- function()			# 判別係数を計算する
	{
		lxi <- lx[1:ni]
		ncase.k <- ncase-ng
		cat("\n***** 判別関数 *****\n")
		for (g1 in 1:(ng-1)) {
			for (g2 in (g1+1):ng) {
				xx <- Mean[lxi, g1]-Mean[lxi, g2]
				fn <- w[lxi, lxi]%*%xx*ncase.k
				fn0 <- -sum(fn*(Mean[lxi, g1]+Mean[lxi, g2])*0.5)
				dist <- sqrt(sum(xx*fn))
				errorp <- pnorm(dist*0.5, lower.tail=FALSE)
				result3 <- data.frame(判別係数= c(fn, fn0), 標準化判別係数=c(sd[lxi]*fn, NA))
				rownames(result3) <- c(vname[lxi], "定数項")
				class(result3) <- c("sdis", "data.frame")	# print.sdis を使うための設定
				cat(sprintf("\n%s と %s の判別\n", grp.lab[g1], grp.lab[g2]))
				cat(sprintf("マハラノビスの汎距離: %.5f\n", dist))
				cat(sprintf("理論的誤判別率:　　　 %s\n", format.pval(errorp, 3, 1e-3)))
				print(result3)
			}
		}
		return(list(fn=fn, fn0=fn0))
	}
	
	proc.predict <- function(ans)		
	{
		nc0 <- 0
		ncase.k <- ncase-ng	
		lxi <- lx[1:ni]
		data <- as.matrix(data)[, lxi, drop=FALSE]	# モデル中の独立変数を順序通りに取り出す
		dis <- matrix(0, ncase, ng)
		for (g in 1:ng) {				# 各群の中心までの距離を計算する
			xx <- t(t(data)-Mean[lxi, g])
			dis[,g] <- rowSums(xx%*%w[lxi, lxi]*xx*ncase.k)
		}
		pred.group <- grp.lab[apply(dis, 1, which.min)]	# 判別された群
		P <- pchisq(dis, p, lower.tail=FALSE)		# その群に属するとしたとき，距離がそれより大きくなる確率
		result4 <- data.frame(実際の群=group, 判別された群 =pred.group,
					正否 =ifelse(group==pred.group, "  ", "##"), dis,
					matrix(format.pval(P, 3, 1e-3), ncase))
		colnames(result4)[4:(3+2*ng)] <- c(paste("二乗距離", 1:ng, sep=""), paste("P値", 1:ng, sep=""))
		if (ng == 2) {					# 判別値を計算するのは2群判別の場合のみ
			fn <- ans$fn				# 判別係数
			fn0 <- ans$fn0				# 定数項
			result4$dfv <- data%*%fn+fn0		# 判別値
			colnames(result4)[8] <- "判別値"
		}
		class(result4) <- c("sdis", "data.frame")	# print.sdis を使うための設定
		result5 <- xtabs(~result4$実際の群+result4$判別された群)	# 判別結果集計表
		temp <- dimnames(result5)
		dimnames(result5) <- list(実際の群=temp[[1]], 判別された群=temp[[2]])
		return(list(result4=result4, result5=result5))
	}

############## 関数本体

	EPSINV <- 1e-6						# 特異行列の判定値

	if (P.out < P.in) {                                     # Pout ≧ Pin でなければならない
		P.out <- P.in
	}
	step <- 0                                               # step.out にて，大域代入される

	p <- ncol(data)						# 説明変数の個数
	if (p == 1) {
		stepwise <- FALSE
	}
	vname <- colnames(data)					# 変数名（なければ作る）
	if (is.null(vname)) {
		vname <- colnames(data) <- paste("x", 1:p, sep="")
	}
	gname <- names(group)
	group <- factor(as.matrix(group))			# 群を表す変数を取り出す（factor にしておく）
	ok <- complete.cases(data, group)			# 欠損値を含まないケース
	data <- as.data.frame(data[ok,])
	group <- group[ok]
	ncase <- nrow(data)					# ケース数
	grp.lab <- levels(group)				# 群の名前
	ng <- nlevels(group)					# 群の個数
	if (ng <= 1) {
		stop("1群しかありません")
	}
	split.data <- split(data, group)
	Mean <- cbind(matrix(sapply(split.data, colMeans), ncol=ng), colMeans(data))
	dimnames(Mean) <- list(vname, c(grp.lab, "全体"))

	num <- c(sapply(split.data, nrow), ncase)

	if (verbose) {
		cat(sprintf("有効ケース数： %i\n", ncase))
		cat(sprintf("群を表す変数： %s\n\n", gname))
		cat("***** 平均値 *****\n")
		print(Mean)
	}
	if (any(num < 2)) {
		stop("ケース数が1以下の群があります")
	}
	t <- var(data)*(ncase-1)
	w <- matrix(colSums(t(matrix(sapply(split.data, var), ncol=ng))*(num[1:ng]-1)), p)
	dimnames(w) <- dimnames(t)
	detw <- dett <- 1

	sd2 <- sqrt(diag(w)/ncase)
	r <- w/outer(sd2, sd2)/ncase
	if (verbose) {
		cat("\n***** プールされた群内相関係数行列 *****\n\n")
		print(r)
	}
	sd <- sqrt(diag(t)/ncase)
	if (stepwise == FALSE) {				# 変数選択をしないとき
		for (ip in 1:p) {
			ans <- sweep.sdis(w, detw)
			w <- ans$r
			detw <- ans$det
			ans <- sweep.sdis(t, dett)
			t <- ans$r
			dett <- ans$det
		}
		lx <- 1:p					# モデルに含まれる説明変数の列番号を保持
		ni <- p						# モデルに含まれる説明変数の個数
		ans.step.out <- step.out(0)
	}
	else {							# 変数選択をするとき
		if (verbose) {
			cat(sprintf("\n変数編入基準    Pin:  %.5g\n",P.in))
			cat(sprintf("変数除去基準    Pout: %.5g\n", P.out))
		}
		lx <- integer(p)				# モデルに含まれる説明変数の列番号を保持
		ni <- 0						# モデルに含まれる説明変数の個数
		while (ni != p) {				# ステップワイズ変数選択
			ans.max <- fmax()			# 編入候補変数を探索
			P <- (ncase-ng-ni)/(ng-1)*ans.max[1]	# F 値から
			P <- pf(P, ng-1, ncase-ng-ni, lower.tail=FALSE) # P 値を求める
			ip <- ans.max[2]			# 変数の位置
			if (verbose) cat(sprintf("編入候補変数: %-15s   P : %s", vname[ip], format.pval(P, 3, 1e-3)))
			if (P > P.in) {
				if (verbose) cat("  ***** 編入されませんでした\n")
				break;				# これ以上の変数は組み込まれない。ステップワイズ選択の終了
			}
			if (verbose) cat("  ***** 編入されました\n")
			ni <- ni+1
			lx[ni] <- ip
			ans <- sweep.sdis(w, detw)
			w <- ans$r
			detw <- ans$det
			ans <- sweep.sdis(t, dett)
			t <- ans$r
			dett <- ans$det
			ans.step.out <- step.out(1)		# 途中結果を出力する

			repeat {				# 変数除去のループ
				ans.min <- fmin()		# 除去候補変数について同じく
				P <- (ncase-ng-ni+1)/(ng-1)*ans.min[1]
				P <- pf(P, ng-1, ncase-ng-ni+1, lower.tail=FALSE)
				ip <- ans.min[2]
				if (verbose) cat(sprintf("\n除去候補変数: %-15s   P : %s", vname[ip], format.pval(P, 3, 1e-3)))
				if (P <= P.out) {
					if (verbose) cat("  ***** 除去されませんでした\n")
					break			# 変数除去の終了
				}
				else {
					if (verbose) cat("  ***** 除去されました\n")
					lx <- lx[-which(lx == ip)]
					ni <- ni-1
					ans <- sweep.sdis(w, detw)
					w <- ans$r
					detw <- ans$det
					ans <- sweep.sdis(t, dett)
					t <- ans$r
					dett <- ans$det
					ans.step.out <- step.out(2) # 途中結果を出力する
				}
			}
		}
	}

	if (ni == 0) {
		warning(paste("条件（ Pin <", P.in, "）を満たす独立変数がありません"))
	}
	else {
		if (verbose) cat("\n===================== 結果 =====================\n")
		cat("\n***** 分類関数 *****\n")
		print(ans.step.out)
		ans.df <- discriminant.function()
		ans.predict <- proc.predict(ans.df)
		if (predict) {
			cat("\n***** 各ケースの判別結果 *****\n\n")
			print(ans.predict$result4)
			cat("\n    メモ:「二乗距離」とは，各群の重心までのマハラノビスの汎距離の二乗です。\n")
			cat("         P値は各群に属する確率です。\n")
		}
		cat("\n***** 判別結果集計表 ****\n\n")
		print(ans.predict$result5)
		ans <- list(分類関数=ans.step.out, 個々の判別=ans.predict$result4, 判別結果集計表=ans.predict$result5)
		class(ans) <- c("sdis", "list")			# plot.sdis を使うための設定
		invisible(ans)
	}
}
# print メソッド
print.sdis <- function(result)					# sdis が返すオブジェクト
{
	if (class(result)[2] == "list") {
		print.default(result)
	}
	else if (class(result)[2] == "data.frame") {
		result <- capture.output(print.data.frame(result, digits=5))
		result <- gsub("$", "\\\n", result)
		result <- gsub("<NA>", "    ", result)
		result <- gsub("NA", "  ", result)
		cat("\n", result, sep="")
	}
}
# plot メソッド
plot.sdis <- function(	result,					# sdis が返すオブジェクト
			which=c("boxplot", "barplot", "scatterplot"),		# 描画するグラフの種類
			nclass=20,				# barplot のときのおよその階級数
			pch=1:2,				# scatterplot を描く記号
			col=1:2,				# scatterplot の記号の色
			xpos="topright", ypos=NULL,			# scatterplot の凡例の位置
			...)					# boxplot, barplot へ引き渡す引数
{
	if (nlevels(result$個々の判別$実際の群) == 2) {
		which <- match.arg(which)
		if (which == "boxplot") {			# boxplot
			plot(result$個々の判別$判別値 ~ result$個々の判別$実際の群, xlab="群", ylab="判別値", ...)
		}
		else if (which == "barplot") { 						# barplot
			tbl <- table(result$個々の判別$実際の群, cut(result$個々の判別$判別値,
					breaks=pretty(result$個々の判別$判別値, n=nclass)))
			barplot(tbl, beside=TRUE, legend=TRUE, xlab="判別値", ...)
		}
		else {						# scatterplot 各群の重心までの二乗距離
			group <- result$個々の判別$実際の群
			group.levels <- levels(group)
			distance1 <- result$個々の判別$二乗距離1
			distance2 <- result$個々の判別$二乗距離2
			max1 <- max(distance1)
			max2 <- max(distance2)
			max0 <- max(max1, max2)
			plot(distance1, distance2, col=col[as.integer(group)], pch=pch[as.integer(group)],
				xlim=c(0, max0), xlab=paste(group.levels[1], "の重心への二乗距離"),
				ylim=c(0, max0), ylab=paste(group.levels[2], "の重心への二乗距離"), asp=1, ...)
			abline(0, 1, lty=2)
			text(max1, max2/2, paste(group.levels[2], "に判別される領域"), pos=2)
			text(0, max2+strheight("H")*1.5, paste(group.levels[1], "に判別される領域"), pos=4)
			legend(x=xpos, y=ypos, legend=group.levels, col=col, pch=pch)
		}
	}
	else {
		warning("3群以上の場合にはグラフ表示は用意されていません")
	}
}
# データ行列から類似度行列を作る
similarity.matrix <- function(	dat,				# データ行列
				power=1,			# 距離のべき乗数（ユークリッド二乗距離なら 2 を指定）
				...)				# dist 関数への引数（距離の種類などの指定）
{
	d <- as.matrix(dist(dat, ...))
	if (power != 1) {
		d <- d^power
	}
	d <- -d
	diag(d) <- 0
	return(d)
}
# シンプレックス法によるパラメータ推定
simplex <- function(	fun,					# 残差平方が最小値となるパラメータを探す目的関数
			start,					# パラメータの初期値ベクトル
			x,					# x 値ベクトル
			y,					# y 値ベクトル
			MAXIT=10000,				# 繰り返し数
			EPSILON=1e-7,				# 推定許容誤差
			LO=0.8, HI=1.2,				# パラメータの初期値ベクトルから 3 組のパラメータベクトルを作るときの倍数
			plot.flag=FALSE,			# TRUE のときには，あてはめ図を描く
			...)					# plot, lines に渡すパラメータ
{
	residual <- function(x, y, p)				# 残差平方和を求める関数
	{
	    	return(sum((y-fun(x, p))^2))
	}
	ip3 <- (ip2 <- (ip1 <- (ip <- length(start))+1)+1)+1
	pa <- matrix(start, nrow=ip, ncol=ip3)
	diag(pa) <- diag(pa)*runif(ip, min=LO, max=HI)
	res <- c(sapply(1:ip1, function(i) residual(x, y, pa[, i])), 0, 0)
	for (loops in 1:MAXIT) {
		res0 <- res[1:ip1]
		mx <- which.max(res0)
		mi <- which.min(res0)
		s <- rowSums(pa[, 1:ip1])
		if (res[mx] < EPSILON || res[mi] < EPSILON || (res[mx]-res[mi])/res[mi] < EPSILON) {
			break
		}
		i <- ip2
		pa[, ip2] <- (2*s-ip2*pa[, mx])/ip
		res[ip2] <- residual(x, y, pa[, ip2])
		if (res[ip2] < res[mi]) {
			pa[, ip3] <- (3*s-(2*ip1+1)*pa[, mx])/ip
			res[ip3] <- residual(x, y, pa[, ip3])
			if (res[ip3] <= res[ip2]) {
				i <- ip3
			}
		}
		else if (res[ip2] > res[mx]) {
			pa[, ip3] <- s/ip1
			res[ip3] <- residual(x, y, pa[, ip3])
			if (res[ip3] >= res[mx]) {
				for (i in which(1:ip1 != mi)) {
					pa[, i] <- (pa[, i]+pa[, mi])*0.5
					res[i] <- residual(x, y, pa[, i])
				}
				i <- 0 # false
			}
			else {
				i <- ip3
			}
		}
		if (i > 0) {
			pa[, mx] <- pa[, i]
			res[mx] <- res[i]
		}
	}
	p <- pa[, mi]
	residuals <- residual(x, y, p)
	if (plot.flag) {
		plot(y ~ x, ...)
		range <- par()$usr
		x <- seq(range[1], range[2], length=1000)
		lines(x, fun(x, p), ...)
	}
	return(list(converge=loops < MAXIT, parameters=p,residuals=residuals))
}
# 因子負荷量の大きさの順に変数を並べ替える
sort.loadings <- function(x)					# factanalが返すオブジェクト
{
	a <- x$loadings
	y <- abs(a)						# 因子負荷量の絶対値
	z <- apply(y, 1, which.max)				# 各変数をどの因子に含めるべきか
	loadings <- NULL					# 結果
	for (i in 1:ncol(y)) {
		b <- a[z == i,, drop=FALSE]
		if (nrow(b)) {
			t <- order(b[, i, drop=FALSE], decreasing=TRUE)	# 因子単位で並べ替え情報を得る
			loadings <- rbind(loadings, b[t,, drop=FALSE])
		}
	}
	class(loadings) <- "loadings"				# クラスの設定
	return(loadings)					# 結果を返す
}
spearman2 <- function(x, y)		# 2変数間のスピアマンの順位相関係数
{
	cor(cbind(x, y), use="complete.obs", method="spearman")[1, 2]
}

spearman <- function(data.matrix)	# スピアマンの順位相関係数行列
{
	cor(data.matrix, use="complete.obs", method="spearman")
}
sreg <- function(	data,					# 最終列が従属変数
			stepwise=TRUE,				# ステップワイズ変数選択をする
			P.in=0.05,				# Pin
			P.out=0.05,				# Pout （Pout ≧ Pin のこと）
			predict=FALSE,				# 各ケースの予測値を出力する
			verbose=TRUE)				# ステップワイズ変数選択の途中結果を出力する
{
	sd2 <- function (x, na.rm=FALSE) 
	{
		if (is.matrix(x)) {
			apply(x, 2, sd, na.rm=na.rm)
		}
    		else if (is.data.frame(x)) {
		        sapply(x, sd, na.rm=na.rm)
		}
	}

	stat <- function()					# 基本統計量の出力
	{
		tbl <- data.frame(平均値=mean, 不偏分散 =var, 標準偏差=sd)
		print(tbl, digits=5)
		if (any(var == 0)) {
			stop("分散が0になる変数があります")
		}
		cat("\n***** 相関係数行列 *****\n")
		print(r, digits=5)
	}

	step.out <- function(isw)				# 変数選択の途中結果を出力
	{
		step <<- step+1					# ステップ数
		syy <- ss[q1]					# 従属変数の平方和
		names(syy) <- NULL				# 後々余計な名前が付くのを防ぐ
		lxi <- lx[1:ni]					# モデルに含まれている変数の列番号
		b1 <- r[q1, lxi]				# 標準化偏回帰係数
		b <- b1*sd[q1]/sd[lxi]				# 偏回帰係数
		b0 <- mean[q1]-sum(b*mean[lxi])			# 定数項
		sr <- sum(b*ss[lxi])				# 回帰の平方和
		se <- syy-sr					# 残差の平方和
		if (se < 0 && abs(se/syy) < 1e-12) {		# 負の値になっても誤差範囲なら許容
			se <- 0
		}
		idf1 <- ni					# 回帰の平方和の自由度
		idf2 <- nc-ni-1					# 誤差の平方和の自由度
		vr <- sr/idf1					# 回帰の平均平方
		ve <- se/idf2					# 残差の平均平方
		if (ve == 0) {					# 完全に当てはまる場合
			f <- fp9 <- NA
			stde <- tv <- pp <- rep(NA, ni)
			seb0 <- tv0 <- pp0 <- NA
			warning("従属変数の予測値と観測値が完全に一致します。分析の指定に間違いはないですか？")
		}
		else {						# 普通はこちら
			f <- vr/ve				# 分散分析の F 値
			fp9 <- pf(f, idf1, idf2, lower.tail=FALSE) # 対応する P 値
		}
		rhoy <- se/syy					# 残差平方和の割合

		if (ve != 0) {
			seb <- r[cbind(lxi, lxi)]*rhoy/idf2*var[q1]/var[lxi]
			stde <- sqrt(seb)			# 偏回帰係数の標準誤差
			tv <- abs(b/stde)			# t 値
			pp <- pt(tv, idf2, lower.tail=FALSE)*2	# P 値
			temp <- mean[lxi]/sd[lxi]
			seb0 <- sqrt((1/nc+temp %*% r[lxi, lxi, drop=FALSE] %*% temp/(nc-1))*ve)
			tv0 <- abs(b0/seb0)
			pp0 <- pt(tv0, idf2, lower.tail=FALSE)*2
		}
		stde[ni+1] <- seb0
		tv[ni+1] <- tv0
		pp[ni+1] <- pp0

		adj.r2 <- 1-ve/syy*(nc-1)			# 自由度調整済み決定係数
		if (rhoy != 0 && isw != 0) {
			r2.change <- old.rhoy-rhoy		# 決定係数の増分
			f.change <- r2.change*idf2/rhoy 	# F 値に換算
			p.change <- if (f.change < 0) NA	# P 値
					else pf(f.change, 1, idf2, lower.tail=FALSE)
		}
		else {
			p.change <- NA
		}
		if (isw != 0 && verbose) {			# ステップのまとめ
			cat(sprintf("\n***** ステップ %i *****   ", step))
			cat(sprintf("%s変数: %s\n", c("編入", "除去")[isw], vname[ip]))
		}
		multico <- r[cbind(lxi, lxi)]			# 分散拡大要因
		result1 <- data.frame(	偏回帰係数=c(b, b0),	# 結果をデータフレームにする
					標準誤差=stde,
					"t値"=tv,
					"P値"=format.pval(pp, 3, 1e-3),
					標準化偏回帰係数=c(b1, NA),
					トレランス=c(1/multico, NA),
					分散拡大要因=c(multico, NA))
		rownames(result1) <- c(vname[lxi], "定数項")	# 行名を付ける
		class(result1) <- c("sreg", "data.frame")	# print.sreg を使うための設定
		if (verbose) print(result1)			# 結果出力
		
		result3 <- data.frame(	平方和=c(sr, se, syy),	# 分散分析表のデータフレーム
					自由度=c(ni, idf2, nc-1),
					平均平方=c(vr, ve, NA),
					"F値"=c(f, NA, NA),
					"P値"=c(format.pval(fp9, 3, 1e-3), NA, NA))
		rownames(result3) <- c("回帰", "残差", "全体")	# 行名を付ける
		class(result3) <- c("sreg", "data.frame")	# print.sreg を使うための設定
		if (verbose) print(result3)			# 結果出力

		loglik <- 0.5*(sum(-nc*(log(2*pi)+1-log(nc)+log(se))))
		AIC <- 2*ni+4-2*loglik
		result4 <- c(重相関係数=sqrt(1-rhoy),
				"決定係数（重相関係数の二乗）"=1-rhoy)
		if (adj.r2 > 0) {
			result4 <- c(result4, 自由度調整済み重相関係数の二乗=adj.r2)
		}
		result4 <- c(result4, 対数尤度=loglik, AIC=AIC)
		class(result4) <- c("sreg", "numeric")		# print.sreg を使うための設定
		if (verbose) print(result4)			# 結果出力

		if (stepwise && !is.na(p.change)) {
			result5 <- c(	決定係数の増分=r2.change,
					"増分に対するF値"=f.change,
					"第1自由度"=1,
					"第2自由度"=idf2,
					"増分に対するP値"=p.change)
			class(result5) <- c("sreg", "numeric")	# print.sreg を使うための設定
			if (verbose) print(result5)		# 結果出力
		}
		else {
			result5 <- NA
		}
		old.rhoy <<- rhoy
		return(list(	b=b, b0=b0, ve=ve, coefficients=result1,
				anova.table=result3, Rs=result4, delta=result5))
	}

	fmax <- function()					# モデルに取り入れる変数の探索
	{
		kouho <- 1:q
		if (ni > 0) {
			kouho <- (1:q)[-lx[1:ni]]
		}
		kouho2 <- cbind(kouho, kouho)
		temp <- 1/(r[kouho2]*r[q1, q1]/(r[q1, kouho]^2)-1)
		ip <- which.max(temp)
		return(c(temp[ip], kouho[ip]))
	}
	
	fmin <- function()					# モデルから捨てる変数の探索
	{
		kouho <- lx[1:ni]
		kouho2 <- cbind(kouho, kouho)
		temp <- r[q1, kouho]^2/(r[kouho2]*r[q1, q1])
		ip <- which.min(temp)
		return(c(temp[ip], lx[ip]))
	}

	sweep.sreg <- function(r)				# 掃き出し法
	{
		ap <- r[ip, ip]
		if (abs(ap) <= EPSINV) {
			stop("正規方程式の係数行列が特異行列です")
		}
		for (i in 1:q1) {
			if (i != ip) {
				temp <- r[ip, i]/ap
				for (j in 1:q1) {
					if (j != ip) {
						r[j, i] <- r[j, i]-r[j, ip]*temp
					}
				}
			}
		}
		r[,ip] <- r[,ip]/ap
		r[ip,] <- -r[ip,]/ap
		r[ip, ip] <- 1/ap
		return(r)
	}

	proc.predict <- function()				# 予測値，標準化残差を求める
	{
		lxi <- lx[1:ni]
		ve <- ans.step.out$ve				# 残差の平均平方
		constant <- ans.step.out$b0			# 定数項
		b <- ans.step.out$b				# 偏回帰係数
		y <- data[,ncol(data)]				# 観察値（従属変数）
		data <- as.matrix(data)[, lxi, drop=FALSE]	# モデル中の独立変数を順序通りに取り出す
		est <- data%*%b+constant			# 予測値
		res <- y-est					# 残差
		r <- r[lxi, lxi]
		data <- scale(data)*sqrt(nc/(nc-1))		# データの標準化
		stres <- apply(data, 1, function(arg) arg%*%r%*%arg)
		stres <- res/sqrt(ve*(1-(stres+1)/nc))		# 標準化残差
		result <- data.frame(観察値=y, 予測値=est, 残差=res, 標準化残差=stres)
		return(result)
	}

############## 関数本体

	EPSINV <- 1e-6						# 特異行列の判定値
	MULTICO <- 10						# 分散拡大要因の基準値
	if (P.out < P.in) {					# Pout ≧ Pin でなければならない
		P.out <- P.in
	}
	step <- 0						# step.out にて，大域代入される
	old.rhoy <- 1						# step.out にて，大域代入される

	vname <- colnames(data)					# 変数名（なければ作る）
	if (is.null(vname)) {
		vname <- colnames(data) <- c(paste("x", 1:(ncol(data)-1), sep=""), "y")
	}
	data <- subset(data, complete.cases(data))		# 欠損値を含むデータを除く
	nc <- nrow(data)					# ケース数
	q1 <- ncol(data)					# 従属変数の位置
	if (verbose) {
		cat(sprintf("有効ケース数： %i\n", nc))
		cat(sprintf("従属変数：　　 %s\n", vname[q1]))
	}

	q <- q1-1						# 独立変数の個数
	if (stepwise == FALSE && nc <= q1) {
		stop("独立変数の個数よりデータ数が多くなければなりません")
	}

	if (q == 1) {						# 単回帰
		stepwise <- FALSE
	}
	mean <- colMeans(data)					# 平均値
	sd <- sd2(data)						# 標準偏差
	var <- sd^2						# 分散
	ss <- (var(data)*(nc-1))[,q1]			# 共変動ベクトル
	r <- cor(data)						# 相関係数行列
	if (verbose) stat()					# 基礎統計量の出力

	stde <- tv <- pp <- numeric(q1)				# メモリ確保
	if (stepwise == FALSE) {				# 変数選択をしないとき
		for (ip in 1:q) {
			r <- sweep.sreg(r)
		}
		lx <- 1:q					# モデルに含まれる独立変数の列番号を保持
		ni <- q						# モデルに含まれる独立変数の個数
		ans.step.out <- step.out(0)
	}
	else {							# 変数選択をするとき
		if (verbose) {
			cat(sprintf("\n変数編入基準    Pin:  %.5g\n", P.in))
			cat(sprintf("変数除去基準    Pout: %.5g\n", P.out))
		}
		lx <- numeric(q)				# モデルに含まれる独立変数の列番号を保持
		ni <- 0						# モデルに含まれる独立変数の個数
		while (ni < min(q, nc-2)) {
			ans.max <- fmax()			# 編入候補変数を探索
			P <- (nc-ni-2)*ans.max[1]		# F 値から
			P <- pf(P, 1, nc-ni-2, lower.tail=FALSE) # P 値を求める
			ip <- ans.max[2]			# 変数の位置
			if (verbose) cat(sprintf("編入候補変数: %-15s   P : %s",vname[ip], format.pval(P, 3, 1e-3)))
			if (P > P.in) {
				if (verbose) cat("  ***** 編入されませんでした\n")
				break
			}
			if (verbose) cat("  ***** 編入されました\n")
			ni <- ni+1
			lx[ni] <- ip
			r <- sweep.sreg(r)
			ans.step.out <- step.out(1)		# 途中結果を出力する

			repeat {				# 変数除去のループ
				ans.min <- fmin()		# 除去候補変数について同じく
				P <- (nc-ni-1)*ans.min[1]
				P <- pf(P, 1, nc-ni-1, lower.tail=FALSE)
				ip <- ans.min[2]
				if (verbose) cat(sprintf("除去候補変数: %-15s   P : %s",vname[ip], format.pval(P, 3, 1e-3)))
				if (P <= P.out) {
					if (verbose) cat("  ***** 除去されませんでした\n")
					break			# 変数除去の終了
				}
				else {
					if (verbose) cat("  ***** 除去されました\n")
					lx <- lx[-which(lx == ip)]
					ni <- ni-1
					r <- sweep.sreg(r)
					ans.step.out <- step.out(2) # 途中結果を出力する
				}
			}
		}
	}
	if (ni == 0) {
		warning(paste("条件（ Pin <", P.in, "）を満たす独立変数がありません"))
	}
	else {
		result6 <- proc.predict() 			# 予測
		ans <- list(	分析結果=ans.step.out$coefficients,
				分散分析表=ans.step.out$anova.table,
				決定係数=ans.step.out$Rs,
				予測=result6)
		if (stepwise) {
			if (verbose) cat("\n===================== 結果 =====================\n")
			print(ans$分析結果)
			print(ans$分散分析表)
			print(ans$決定係数)
			if (predict) {
				cat("\n")
				print(ans$予測, digits=5)
			}
		}
		class(ans) <- c("sreg", "list")			# plot.sdis を使うための設定
		invisible(ans)
	}
}

print.sreg <- function(result)					# sreg クラスのオブジェクトの print メソッド
{
	if (class(result)[2] == "list") {
		print.default(result)
	}
	else if (class(result)[2] == "data.frame") {
		result <- capture.output(print.data.frame(result, digits=5))
		result <- gsub("$", "\\\n", result)
		result <- gsub("<NA>", "    ", result)
		result <- gsub("NA", "  ", result)
		cat("\n", result, sep="")
	}
	else if (length(result) == 5) {
		nm <- names(result)
		if (nm[1] == "重相関係数") {
			cat(sprintf(	"\n%s\t%.5f\n%s\t%.5f\n%s\t%.5f\n%s\t%.5f\n%s  \t%.5f\n",
					nm[1], result[1], nm[2], result[2], nm[3], result[3],
					nm[4], result[4], nm[5], result[5]))
		}
		else {
			cat(sprintf(	"\n%s\t%.5f\n%s\t%.5f\n%s\t%i\n%s\t%i\n%s\t%.5f\n\n",
					nm[1], result[1], nm[2], result[2], nm[3], result[3],
					nm[4], result[4], nm[5], result[5]))
		}
	}
}

plot.sreg <- function(result, which=c("stdres", "qqplot", "slope1"), ...) # sdis クラスのオブジェクトの plot メソッド
{
	which <- match.arg(which)
	if (which == "stdres") {				# 予測値--標準化残差
		plot(result$予測$標準化残差 ~ result$予測$予測値, xlab="予測値", ylab="標準化残差", ...)
	}
	else if (which == "qqplot") {				# 標準化残差の Q-Q プロット
		n <- length(result$予測$標準化残差)
		qqnorm(result$予測$標準化残差, xlab="理論値", ylab="標準化残差", ...)
		qqline(result$予測$標準化残差, lty=3)
	}
	else {							# 予測値--観察値
		plot(result$予測$観察値 ~ result$予測$予測値, xlab="予測値", ylab="観察値", asp=1, ...)
		abline(0, 1, lty=3)
	}
}

# 定点を通る直線回帰式の傾き
sregc <- function(	x,					# 独立変数			
			y,					# 従属変数
			cxy=NULL)				# 定点の y 座標
{
	OK <- complete.cases(x, y)				# 欠損値を持つケースを除く
	x <- x[OK]
	y <- y[OK]
	n <- length(x)						# サンプルサイズ
	sxy <- sum(x*y)						
	sxx <- sum(x^2)
	sx <- sum(x)
	sy <- sum(y)
	if (is.null(cxy)) {					# cxy が NULL なら，平均値
		cx <- mean(x)
		cy <- mean(y)
	}
	else {							# 2 つのベクトルとして指定されていれば
		cx <- cxy[1]
		cy <- cxy[2]
	}
	return((sxy-cy*sx-cx*sy+n*cx*cy)/(sxx-2*cx*sx+n*cx^2))
}
# 重回帰分析の標準化残差
std.res <- function(	x,					# 説明変数だけのデータ行列
			y)					# 従属変数のデータベクトル
{
	ans <- lm(y ~ x)					# 重回帰分析を行う
	x <- cbind(1, as.matrix(x))
        s2 <- sum(ans$residuals^2)/(nrow(x)-ncol(x))
        h <- sqrt(s2*(1-rowSums((x %*% solve(t(x)%*%x))* x)))
	retv <- cbind(y, ans$fitted.values, ans$residuals, ans$residuals/h)
	colnames(retv) <- c("y","fitted","residual","std.residual")
	return(retv)
}
# dual クラス のための summary メソッド　（dual, pc.dual, ro.dual が利用する）
summary.dual <- function(	x,				# dual が返すオブジェクト
				nf=ncol(x[[1]]),		# 出力する解の数
				weighted=FALSE,			# 相関比で重み付けした解を出力するなら TRUE
				digits=3)			# 出力する数値の小数点以下の桁数
{
	suf <- if (weighted) 4 else 2				# 相関比で重み付けした解も選べる
	str <- if (weighted) "weighted " else ""
	print(round(x[[1]][, 1:nf, drop=FALSE], digits=digits))
	cat(sprintf("\n%srow score\n", str))
	print(round(x[[suf]][, 1:nf, drop=FALSE], digits=digits))
	cat(sprintf("\n%scolumn score\n", str))
	print(round(x[[suf+1]][, 1:nf, drop=FALSE], digits=digits))
}
# 分割表形式で与えられたデータに基づいて，ケンドールのτb（ケンドールの順位相関係数）を計算する
tau.b <- function(f)							# 分割表（合計欄を含めない）
{
	R <- nrow(f)							# 行数
	C <- ncol(f)							# 列数
	n <- sum(f)							# 全数
	rt <- rowSums(f)						# 行和
	ct <- colSums(f)						# 列和
	dr <- n^2-sum(rt^2)
	dc <- n^2-sum(ct^2)
	g <- matrix(0, nr=R+2, nc=C+2)
	cada <- f
	g[2:(R+1), 2:(C+1)] <- f
	PQ <- 0
	for (i in 1:R) {
		for (j in 1:C) {
			cada[i, j] <- 	sum(g[1:i, 1:j], g[(i+2):(R+2), (j+2):(C+2)]) -
					sum(g[1:i, (j+2):(C+2)], g[(i+2):(R+2), 1:j])
			PQ <- PQ+g[i+1, j+1]*cada[i, j]
		}
	}
	taub <- PQ/sqrt(dr*dc)						# τb
	ase0 <- 2*sqrt((sum(f*cada^2)-PQ^2/n)/(dr*dc))			# 標準誤差
	Z <- taub/ase0							# 検定統計量
	P <- pnorm(abs(Z), lower.tail=FALSE)*2				# P 値
	for (i in 1:R) {
		for (j in 1:C) {
			f[i, j] <- f[i, j]*(2*sqrt(dr*dc)*cada[i, j]+taub*(rt[i]*dc+ct[j]*dr))^2
		}
	}
	ase1 <- sqrt(sum(f)-n^3*taub^2*(dr+dc)^2)/(dr*dc)		# 別法（よく使われる方法）
	Z2 <- taub/sqrt((4*n+10)/(9*n*(n-1)))				# ケンドールの順位相関係数の検定のときに出てくる式
	P2 <- pnorm(abs(Z2), lower.tail=FALSE)*2			# P 値
	c("tau b"=taub, "ase1"=ase1, "ase0"=ase0, "Z value"=Z, "P value"=P, "Z value 2"=Z2, "P value 2"=P2)
}
# 分割表データから，元のデータを構成し，R にある cor.test を使う。
tau.b2 <- function(f)
{
	x <- rep(row(f), f)
	y <- rep(col(f), f)
	cor.test(x, y, method="kendall")
}
# 分割表を与えて，その分割表が得られる元のデータ（二変数データ）を再現する
tenkai <- function(f)
{
	list(x=rep(row(f), f), y=rep(col(f), f))
}
# 分割表を与えて，その分割表が得られる元のデータ（二変数データ）を再現する
tenkai <- function(f)
{
	list(x=rep(row(f), f), y=rep(col(f), f))
}
# テューキーの方法による線形比較
tlc <- function(n,				# 各群のデータ個数のベクトル
		m,				# 各群の平均値のベクトル
		u,				# 各群の不偏分散のベクトル
		g1,				# 第 1 グループの指定
		g2,				# 第 2 グループの指定
		conf.level=0.95)		# 線形比較の信頼区間を計算する信頼率
{
	stopifnot(	length(n) == length(m),	# データのチェック
			length(m) == length(u),	# n, m, u の要素数は同じでなくてはならない
			n > 1,			# 各群の例数は 2 以上でなくてはならない
			u > 0,			# 不偏分散は正の値でなくてはならない
			floor(n) == n,		# 例数は整数でなくてはならない
			n == n[1],		# 各群の例数は同じでなくてはならない
			floor(g1) == g1,	# 群の番号は整数でなくてはならない
			floor(g2) == g2)	# 群の番号は整数でなくてはならない
	method <- "テューキーの方法による線形比較"
	data.name <- paste(deparse(substitute(g1)), "and", deparse(substitute(g2)))
	ng <- length(n)				# 群の数
	nc <- sum(n)				# 全例数
	dfw <- nc-ng				# 群内分散の自由度
	Vw <- sum(u*(n-1))/dfw			# 群内分散
	n1 <- length(g1)			# 第 1 グループに含まれる群の数
	n2 <- length(g2)			# 第 2 グループに含まれる群の数
	g0 <- (1:ng)[-c(g1, g2)]		# どちらのグループにも含まれない群の番号
	n0 <- ng-n1-n2				# どちらのグループにも含まれない群の数
	weight <- rep(c(1/n1, -1/n2, 0),	# 線形比較の重み
			c(n1, n2, n0))[order(c(g1, g2, g0))]
	theta <- sum(weight*m)			# 線形比較 θ
	sq <- sqrt(Vw/n[1])			# θの標準誤差
	q <- qtukey(1-conf.level, ng, dfw,	# αに対応するステューデント化した範囲
			lower.tail=FALSE)
	p <- ptukey(abs(theta/sq), ng, dfw,	# P 値
			lower.tail=FALSE)
	conf.int <- theta-c(1, -1)*q*sq		# θの信頼区間
	return(structure(list(statistic=c(theta=theta), parameter=c(a=ng, df=dfw), p.value=p,
		conf.int=conf.int, method=method, data.name=data.name, contrast=list(g1, g2)),
		class="htest"))
}
# 重回帰分析の回帰診断の一つとして，トレランスを計算する
tolerance <- function(x)					# 説明変数だけのデータ行列
{
	if (is.null(colnames(x))) {			        # 名前が付いていないときには仮の名前を付ける
		colnames(x) <- paste("Var", 1:ncol(x), sep="")
	}
	x <- subset(x, complete.cases(x))			# 欠損値を持つケースを除く
	VIF <- diag(solve(cor(x)))				# 分散拡大要因（相関係数行列の逆行列の対角成分）
	tolerance <- 1/VIF					# トレランス（VIF の逆数）
	result <- data.frame(tolerance, VIF)			# 結果をデータフレームにする
	return(result)
}
# 三角行列の要素を与えて対称行列を作る
tri.mat <- function(	x,					# 対角成分を含む三角行列の要素
			n=floor(sqrt(1+8*length(x))-1)/2,	# 指定不要（対称行列のサイズ）
			byrow=TRUE,				# x が行優先で使用されるとき TRUE
			lower=TRUE)				# x が下三角行列のとき TRUE
{
	stopifnot(length(x) == n*(n+1)/2)			# 要素数に過不足がないかチェック
	mat <- diag(n)						# 正方行列を作る
	mat[(if (xor(byrow, lower)) lower.tri else upper.tri)(mat, diag=TRUE)] <- x
	mat <- t(mat)+mat
	diag(mat) <- diag(mat)/2
	return(mat)
}
# Tukey の方法による多重比較
# Games-Howell の方法も選択できるように拡張 2009/08/03
tukey <- function(	data,					# 観察値ベクトル
			group,					# 群変数ベクトル
			method=c("Tukey", "Games-Howell"))	# 手法の選択
{
	OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
	data <- data[OK]
	group <- factor(group[OK])
	n <- tapply(data, group, length)			# 各群のケース数
	a <- length(n)						# 群の数
	phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
	Mean <- tapply(data, group, mean)			# 各群の平均値
	Variance <- tapply(data, group, var)			# 各群の不偏分散
	result1 <- cbind(n, Mean, Variance)			# 各群の統計量
	rownames(result1) <- paste("Group", 1:a, sep="")
	method <- match.arg(method)
	if (method == "Tukey") {				# Tukey の方法
		v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
		t <- combn(a, 2, function(ij)			# 対比較
					abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
		p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
		Tukey <- cbind(t, p)					# 対比較の結果
		rownames(Tukey) <- combn(a, 2, paste, collapse=":")
		return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
	}
	else {							# Games-Howell の方法
		t.df <- combn(a, 2, function(ij) {		# 対比較
					t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
					df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
					return(c(t, df))} )
		t <- t.df[1,]
		df <- t.df[2,]
		p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	# 有意確率を計算する
		Games.Howell <- cbind(t, df, p)			# 対比較の結果
		rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
		return(list(result1=result1, Games.Howell=Games.Howell))
	}
}
#####
#
# 独立 k 標本の平均値，標準偏差を求め，必要なら平均値・代表値の差の検定を行う
#
#####

indep.sample <- function(i,							# 分析対象の変数が入っているデータフレーム上の列番号または変数名ベクトル
			 j,							# 群を表す変数が入っているデータフレーム上の列番号または変数名ベクトル
			 df,							# データフレーム
			 latex=TRUE,						# LaTeX 形式で集計表を出力する（デフォルトは LaTeX 形式）
                         captions=NULL,						# latex=TRUE のときに，各表の表題を表す文字列ベクトルを指定できる（NULL のときはデフォルトの表題）
                         labels=NULL,						# latex=TRUE のときに，各表の label を表す文字列ベクトルを指定できる（NULL のときは付けない）
			 test=c("none", "parametric", "non-parametric"),	# デフォルト none では検定を行わない。検定を行うときはその種類を指定する
			 statistics=c("mean", "median"),			# （平均値，標準偏差）を求めるか（中央値，四分偏差）を求めるかを指定する
			 var.equal=FALSE,					# parametric の場合に等分散性を仮定するかどうかの引数
			 digits=3,						# 平均値，標準偏差を表示するときの小数点以下の桁数
			 output="",						# ファイルに出力するときはファイル名（デフォルトはコンソールに出力）
			 encoding=getOption("encoding"))			# ファイルに出力するときのエンコーディング（デフォルトは OS による）
{

# 下請け関数

	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}

	SIQ <- function(x) return(diff(fivenum(x)[c(2,4)]))                     # 四分偏差を求める下請け関数

	indep.sample.sub <- function(ii, jj)
	{
		group <- colnames(df)[jj]					# 群を表す変数の名前
		df2 <- df[, c(ii, jj)]						# データフレームの列番号 ii, jj から 2 変数を取り出す
		df2 <- subset(df2, complete.cases(df2))				# 欠損値を持つケースを除く
		x <- df2[, 1]							# 分析対象変数
		g <- df2[, 2]							# 群変数
		lst <- list(g)							# by 関数を適用するために，群をリスト化する
		nt <- length(x)							# 全体のサンプルサイズ
		mt <- MEAN(x)							# 全体の平均値
		st <- SD(x)							# 全体の標準偏差
		n <- by(x, lst, length)						# 各群のサンプルサイズ
		m <- by(x, lst, MEAN)						# 各群の平均値（中央値）
		s <- by(x, lst, SD)						# 各群の標準偏差（四分偏差）
		nr <- length(table(lst))					# 群の数
		if (latex) {							# LaTeX 形式で集計結果を出力する
			cat("\n\\begin{table}[htbp]\n", file=output)		# \begin{table}[htbp]
			if (is.null(captions)) {
				cat(sprintf("\\caption{%s別の%sの集計}\n",	# \caption{○○別の□□の集計}
				    group, colnames(df2)[1]), file=output)
			}
			else {
				cat(sprintf("\\caption{%s}\n", captions[index]), file=output)	# \caption{○○○○}
			}
			if (!is.null(labels)) {
				cat(sprintf("\\label{%s}\n", labels[index]), file=output)	# \labels{○○○○}
			}
			cat("\\centering\n", file=output)			# \centering
			cat("\\begin{tabular}{lccc} \\hline\n",file=output)	# \begin{tabular}{lccc} \hline
			cat(sprintf("& \\multicolumn{3}{c}{%s}\\\\ \\cline{2-4}\n", # \multicolumn{3}{c}{□□} \\ \cline{2-4}
			    colnames(df2)[1]), file=output)
			cat(group, "データ数", M.str, S.str, sep=" & ",		# ○○ & データ数 & 平均値 & 標準偏差
			    file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			for (l in 1:nr) {					# 各群について，
				cat(names(n)[l], n[l], sprintf(format, m[l]),	# 群の名前 & 数値 & 数値 & 数値
				    sprintf(format, s[l]), sep=" & ", file=output)
				cat("\\\\", file=output)			# \\
				if (l == nr) cat("\\hline\n", file=output)	# 最後の群なら \hline そうでなければ何もなし
				else cat("\n", file=output)
			}
			cat("全体", nt, sprintf(format, mt),			# 全体 & 数値 & 数値 & 数値
			    sprintf(format, st), sep=" & ", file=output)
			cat("\\\\ \\hline\n", file=output)			# \\ \hline
			cat("\\end{tabular}\n", file=output)			# \end{tabular}
		}
		else {								# tab で区切って出力する
			cat("\n表　", group, "別の", colnames(df2)[1],		# 表　○○別の□□の集計
			    "の集計", sep="", file=output)
			cat("\n", colnames(df2)[1], sep="\t", file=output,	# 　　□□
			    fill=TRUE)
			cat(group, "データ数", M.str, S.str, sep="\t",		# ○○　データ数　平均値　標準偏差
			    file=output, fill=TRUE)
			for (l in 1:nr) {					# 各群について，
				cat(names(n)[l], n[l], sprintf(format, m[l]),	# 群の名前　数値　数値　数値
				    sprintf(format, s[l]), sep="\t",
				    file=output, fill=TRUE)
			}
			cat("全体", nt, sprintf(format, mt),			# 全体　数値　数値　数値
			    sprintf(format, st), sep="\t",
			    file=output, fill=TRUE)
		}
		if (nr == 2) {							# 2 群の場合には，
			if (latex && test != "none") {				# LaTeX 形式で出力するときには，改行して次の行に出力準備
				cat("\\\\ \\noindent\n", file=output)
			}
			if (test == "parametric") {				# 平均値の差の検定のために t.test 関数を使う
				res <- t.test(x~g, var.equal=var.equal)		# t.test を呼ぶ
				cat(sprintf(if (latex) "$t$値 = %.3f, 自由度 = %.3f, $P$値 = %.3f\n"
					    else "t 値 = %.3f, 自由度 = %.3f, P 値 = %.3f\n",
					    res$statistic, res$parameter, res$p.value), file=output)
			}
			else if (test == "non-parametric") {			# マン・ホイットニーの U 検定
				res <- wilcox.test(x~g)				# wilcox.test を呼ぶ
				cat(sprintf(if (latex) "$U$ = %.3f, $P$値 = %.3f\n"
					    else "U = %.3f, P 値 = %.3f\n",
					    res$statistic, res$p.value), file=output)
			}
		}
		else if (nr >= 3) {
			if (latex && test != "none") {				# LaTeX 形式で出力するときには，改行して次の行に出力準備
				cat("\\\\ \\noindent\n", file=output)
			}
			if (test == "parametric") {				# 一元配置分散分析
				res <- oneway.test(x ~ g, var.equal=var.equal)
				cat(sprintf(if (latex) "$F$値 = %.3f, 自由度 = (%i, %.3f), $P$値 = %.3f\n"
					    else "F 値 = %.3f, 自由度 = (%i, %.3f), P 値 = %.3f\n",
					    res$statistic, res$parameter[1], res$parameter[2], res$p.value), file=output)
			}
			else if (test == "non-parametric") {			# クラスカル・ウォリス検定
				res <- kruskal.test(x~g)	
				cat(sprintf(if (latex) "$\\chi_{kw}^2$ = %.3f, 自由度 = %i, $P$値 = %.3f\n"
					    else "カイ二乗値(kw) = %.3f, 自由度 = %i, P 値 = %.3f\n", 
					    res$statistic, res$parameter, res$p.value), file=output)
			} 
		}
		if (latex) {							# LaTeX 形式で集計結果を出力する場合は，
			cat("\\end{table}\n", file=output)			# \end{table}
		}
	}

# 関数本体
	if (output != "") {							# 結果をファイルに出力する場合の処理
		output <- file(output, open="w", encoding=encoding)
	}

	test <- match.arg(test)							# 引数が省略形で与えられたときに，正確な名前をとる
	statistics <- match.arg(statistics)					# 引数が省略形で与えられたときに，正確な名前をとる
	if (statistics == "mean") {
		MEAN <- mean							# 位置の母数を求める関数：平均値
		SD <- sd							# 散布度を求める関数：標準偏差
		M.str <- "平均値"
		S.str <- "標準偏差"
	}
	else {
		MEAN <- median							# 位置の母数を求める関数：中央値
		SD <-  SIQ							# 散布度を求める関数：四分偏差
		M.str <- "中央値"
		S.str <- "四分偏差"
	}
	format <- paste("%.", digits, "f", sep="")				# 小数点以下 3 桁で出力する書式
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}

	index <- 0
	for (jj in j) {								# j はベクトルまたはスカラーで，群を表す変数がある列番号
		for (ii in i) {							# i はベクトルまたはスカラーで，分析対象変数がある列番号
			if (ii != jj) {						# i, j の全ての組み合わせについて（ii と jj が違うときのみ），
				index <- index+1				# 集計のための下請け関数 indep.sample.sub を呼ぶ
				indep.sample.sub(ii, jj)
			}
		}
	}

	if (output != "") {							# 結果をファイルに出力した場合の後始末
		close(output)
	}
}
#####
#
# 分析対象変数が，共に数値変数である場合には散布図，何れかが factor である場合には箱髭図を描く
#
#####

twodim.plot <- function(i,							# x 軸に取る変数が入っているデータフレーム上の列番号または変数名ベクトル
			j,							# y 軸に取る変数が入っているデータフレーム上の列番号または変数名ベクトル
			df,							# データフレーム
			k=NULL,							# もしあれば，群を表す変数（factor）が入っているデータフレーム上の列番号または変数名ベクトル
			lm=FALSE,						# 回帰直線を描き込むとき TRUE にする（デフォルト FALSE では，描き入れない）
			cor=c("none", "pearson", "kendall", "spearman"),	# 計算する相関係数を指定する（デフォルト none では，計算しない）
			digits=3,						# 相関係数の小数点以下の桁数
			plot="",						# 散布図を描き出すファイル名（デフォルトは Quarts デバイスに出力）
			type=c("pdf", "png", "jpeg", "bmp", "tiff"),		# 画像フォーマット（plot と併せてファイル名の拡張子として使う）
			width=500,						# 画像の横幅のピクセル数（デフォルトは500ピクセル）
			height=375,						# 画像の高さのピクセル数（デフォルトは375ピクセル）
			xlab=NULL,						# x 軸のラベル（デフォルトは対象変数名）。何も描かないときには空文字列を指定する
			ylab=NULL,						# y 軸のラベル（デフォルトは対象変数名）。何も描かないときには空文字列を指定する
			...)							# 作図関数に渡されるその他の引数
{
	twodim.plot2 <- function(i, j, k)					# 下請け関数。i, j はスカラー
	{
		if (is.null(k)) {
			df2 <- df[, c(i, j)]					# データフレームの列番号 i, j から 2 変数を取り出す
		}
		else {
			if (is.character(k[1])) {
				k <- getNum(k, df)
			}
			df2 <- df[, c(i, j, k )]
		}
		df2 <- subset(df2, complete.cases(df2))				# 欠損値を持つケースを除く
		xlab2 <- if (is.null(xlab)) colnames(df2)[1] else xlab		# x 軸のラベルが指定されていないときには変数名を使う
		ylab2 <- if (is.null(ylab)) colnames(df2)[2] else ylab		# y 軸のラベルが指定されていないときには変数名を使う
		if (is.factor(df2[,1]) && is.factor(df2[,2])) {			# 2 変数ともに数値変数でない場合にはエラー
			cat(i, "列と", j, "列の変数は，共に数値変数ではありません\n")
		}
		else if (is.numeric(df2[,1]) && is.numeric(df2[,2])) {		# 2 変数ともに数値変数の場合には，散布図（および回帰直線など）を描く
			if (is.null(k)) {
				pch <- 1
			}
			else {
				pch <- as.integer(df2[,3])
			}
			plot(df2[,1], df2[,2], xlab=xlab2, ylab=ylab2,		# まずは散布図を描く
			     pch=pch, ...)
			method.name <- switch(cor,				# 引数 cor と変数 method.name の対応付け
						none = "",
						pearson = "ピアソンの積率相関係数",
						kendall = "ケンドールの順位相関係数",
						spearman = "スピアマンの順位相関係数")
			if (method.name == "") {				# 相関係数が不要なら，
				r <- ""						# r は空
			}
			else {							# 3種の相関係数の何れかを計算し，r に文字列として設定する
				r <- paste(method.name, "=", sprintf(format, cor(df2, use="pairwise.complete.obs", method=cor)[1,2]), sep="")
			}
			if (lm) {						# 回帰直線を描き込むなら，
				ans <- lm(df2[,2]~df2[,1])			# 直線回帰分析を行い，
				abline(ans)					# 散布図に回帰直線を描き込み，切片・傾きを str に文字列として設定する
				str <- paste("切片=", sprintf("%g", ans$coefficients[1]), "　傾き=", sprintf("%g", ans$coefficients[2]), sep="")
			}
			else {
				str <- ""					# 回帰直線を描き込まないなら，str は空
			}
			if (r != "" || str != "") {				# 相関係数，または，切片・傾きを書き込むときには，
				str <- paste(r, str, sep="　")			# 両方の結果を str に設定する
				rangex <- range(df2[,1])			# x 軸に取った変数の範囲を求める
				rangey <- range(df2[,2])			# y 軸に取った変数の範囲を求める
				old <- par(xpd=TRUE)				# プロット領域の外にも描き込めるように設定
				text(mean(rangex), 1.1*rangey[2]-0.1*rangey[1], # 適当な位置に str を描く
				     str, pos=3, ...)
				par(old)					# 元に戻す
			}
		}
		else if (is.factor(df2[,1]) && is.numeric(df2[,2])) {		# i 列の変数が factor で，j 列の変数が数値変数なら，
			boxplot(df2[,2]~df2[,1], xlab=xlab2, ylab=ylab2, ...)	# 垂直な boxplot を描く
		}
		else if (is.numeric(df2[,1]) && is.factor(df2[,2])) {		# i 列の変数が数値変数で，j 列の変数が factor なら，
			boxplot(df2[,1]~df2[,2], xlab=xlab2, ylab=ylab2,	# 水平な boxplot を描く
				horizontal=TRUE, ...)
		}
	}
	getNum <- function(str, df) {						# 変数名から列番号を得る
		names <- colnames(df)
		seq_along(names)[names %in% str]
	}
# twodim.plot 関数本体
	cor <- match.arg(cor)							# 引数で短縮形で指定された場合にも cor を正式なものに設定する
	format <- paste("%.", digits, "f", sep="")				# 相関係数の小数点以下の桁数を設定する
	if (plot != "") {							# グラフをファイルに出力するとき plot はファイル名（拡張子を除く）
		type <- match.arg(type)						# 画像ファイルの形式
		if (type == "pdf") {						# pdf なら，一つ一つの画像を別々のファイルに出力するために onefile = FALSE にする
			pdf(sprintf("%s%%03i.pdf", plot), onefile=FALSE,	# pdf は，画像の大きさの指定がインチ単位なので 72dot/inch で換算
			    width=width/72, height=height/72)
		}
		else if (type == "png") {
			png(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else if (type == "bmp") {
			bmp(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else if (type == "jpeg") {
			jpeg(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
		else { # type == "tiff"
			tiff(sprintf("%s%%03i.%s", plot, type), width=width, height=height)
		}
	}
	if (is.character(i[1])) {
		i <- getNum(i, df)
	}
	if (is.character(j[1])) {
		j <- getNum(j, df)
	}
	for (ii in i) {
		for (jj in j) {
			if (ii != jj) {
				twodim.plot2(ii, jj, k)				# 第 1，第 2 引数の全ての組み合わせで図を描くために twodim.plot2 を呼ぶ
			}
		}
	}
	if (plot != "") {							# ファイルに出力しているなら，
		dev.off()							# デバイスを閉じる
	}
}
# 二元配置分散分析
twoway.anova <- function(	x,			# データベクトル
				a,			# 要因 A の factor ベクトル
				b)			# 要因 B の factor ベクトル
{
	name.a <- deparse(substitute(a))
	name.b <- deparse(substitute(b))
	name.ab <- paste(name.a, name.b, sep="*")
	OK <- complete.cases(x, a, b)
	x <- x[OK]
	a <- a[OK]
	b <- b[OK]
	tbl <- table(a, b)				# 各水準の組合せの繰り返し数
	n <- length(x)					# 全データ数
	dft <- n-1					# 自由度
	na <- nrow(tbl)					# 要因 A の水準数
	dfa <- na-1					# 自由度
	nb <- ncol(tbl)					# 要因 B の水準数
	dfb <- nb-1					# 自由度
	gm <- mean(x)					# 全データの平均値
	MSt <- var(x)					# 不偏分散
	sst <- dft*MSt					# 全変動は不偏分散に自由度を掛けたもの
	ma <- tapply(x, a, mean)			# 要因 A の水準ごとの平均値
	mb <- tapply(x, b, mean)			# 要因 B の水準ごとの平均値

	hirei <- function(tbl)				# 繰り返し数が周辺度数に比例するなら TRUE を返す関数
	{
		for (i in 2:nrow(tbl)) {
			temp <- tbl[i,]/tbl[1,]		# 繰り返し数のクロス集計表の，2 行目以降について，1行目のベクトルとの比を計算
			if (any(temp != temp[1])) {	# どれか一つでも比が違えば，比例していない
				return(FALSE)
			}
		}
		return(TRUE)
	}

	if (any(tbl == 0)) {
		result <- "いくつかの水準の組合せにおいて，データがないセルがあります"
	}
	else if (all(tbl == 1)) {			# 繰り返し数が１の場合
		ssa <- nb*sum((ma-gm)^2)		# 要因 A による変動
		ssb <- na*sum((mb-gm)^2)		# 要因 B による変動
		sse <- sst-ssa-ssb			# 誤差変動
		dfe <- dft-dfa-dfb			# 自由度
		ss <- c(ssa, ssb, sse, sst)		# 平方和
		df <- c(dfa, dfb, dfe, dft)		# 自由度
		ms <- ss/df				# 平均平方
		f <- ms/ms[3]				# F 値
		p <- pf(f, df, df[3], lower.tail=FALSE)	# P 値
		f[3:4] <- p[3:4] <- NA			# 計算しないセル
		result <- cbind(ss, df, ms, f, p)
		colnames(result) <- c("SS", "d.f.", "MS", "F value", "P value")
		rownames(result) <- c(name.a, name.b, "e", "T")
	}
	else if(all(tbl >= 2) && hirei(tbl)) {		# 繰り返し数が周辺度数に比例する場合
		mab <- matrix(0, na, nb)		# 水準の組み合わせごとの平均値を計算する
		for (i in 1:na) {
			for (j in 1:nb) {
				mab[i, j] <- mean(x[a == i & b == j])
			}
		}
		sse <- 0				# 残差平方和
		for (i in 1:na) {
			for (j in 1:nb) {
				sse <- sse+sum((x[a == i & b == j]-mab[i,j])^2)
			}
		}
	
		ssa <- sum(rowSums(tbl)*(ma-gm)^2)	# 要因 A による平方和
		ssb <- sum(colSums(tbl)*(mb-gm)^2)	# 要因 B による平方和
		ss <- c(ssa, ssb, sst-ssa-ssb-sse, sse, sst)
		dfab <- dfa*dfb
		df <- c(dfa, dfb, dfab, dft-dfa-dfb-dfab, dft)
		ms <- ss/df				# 平均平方
		f <- ms/ms[4]				# F 値
		p <- pf(f, df, df[4], lower.tail=FALSE)	# P 値
		f[4:5] <- p[4:5] <- NA			# 計算の不要なセル
		result <- cbind(ss, df, ms, f, p)	# 母数モデル
		f[1:2] <- ms[1:2]/ms[3]			# F 値
		p[1:2] <- pf(f[1:2], df[1:2], df[3], lower.tail=FALSE)	# P 値
		result2 <- cbind(ss, df, ms, f, p)	# 変量モデル
		colnames(result) <- colnames(result2) <- c("SS", "d.f.", "MS", "F value", "P value")
		rownames(result) <- rownames(result2) <- c(name.a, name.b, name.ab, "e", "T")
		result <- list(Model1=result, Model2=result2)
	}
	else {						# 一般の場合
	
		mreg <- function(dat)			# 回帰分析を行い，回帰変動，誤差変動，全変動を返す関数
		{
			nc <- ncol(dat)
			ans <- lm(dat[,nc] ~ dat[, -nc])
			St <- var(dat[, nc])*(nrow(dat)-1)
			Se <- sum(ans$residuals^2)
			return(c(St-Se, Se, St))	# 回帰変動，誤差変動，全変動
		}
	
		d1 <- model.matrix(~factor(a))[,-1]	# 要因 A に関するダミー変数
		d2 <- model.matrix(~factor(b))[,-1]	# 要因 B に関するダミー変数
		d3 <- model.matrix(~factor(a)*factor(b))[,-1]	# 要因 A，B の交互作用も含むダミー変数

		r.a.b.ab <- mreg(cbind(d3, x))		# フルモデル
		r.a.b <- mreg(cbind(d1, d2, x))		# 要因 A，B の主効果のみを含むモデル
		r.a <- mreg(cbind(d1, x))		# 要因 A を含むモデル
		r.b <- mreg(cbind(d2, x))		# 要因 B を含むモデル

		ss <- c(r.a.b[1]-r.b[1], r.a.b[1]-r.a[1], r.a.b.ab[1]-r.a.b[1], r.a.b.ab[2], r.a.b.ab[3])	# 変動の構成
		df <- c(dfa, dfb, dfa*dfb, n-na*nb, n-1)
		ms <- ss/df				# 平均平方
		f <- ms/ms[4]				# F 値
		p <- pf(f, df, df[4], lower.tail=FALSE)	# P 値
		f[4:5] <- p[4:5] <- NA			# 計算不要のセル
		result <- cbind(ss, df, ms, f, p)
		colnames(result) <- c("SS", "d.f.", "MS", "F value", "P value")
		rownames(result) <- c(name.a, name.b, name.ab, "e", "T")
	}
	return(result)
}
# 不偏標準偏差（不偏分散の平方根ではない）
unbiased.sd <- function(x)
{
	x <- x[!is.na(x)]					# 欠損値データを除く
	n <- length(x)						# 有効なデータの個数
	return(sd(x)*sqrt((n-1)/2)*gamma((n-1)/2)/gamma(n/2))
}
# 馬鹿馬鹿しいが，Excel にある関数を R で定義してみる
avedev <- function(x)
{
	x <- x[!is.na(x)]
	mean(abs(x-mean(x)))				# 算術平均値からの偏差の絶対値の算術平均値
}

average <- mean	# 単に名前の違い

count <- length	# 単に名前の違い

devsq <- function(x)
{
	x <- x[!is.na(x)]
	(length(x)-1)*var(x)				# 平方和は不偏分散から元に戻す
}

geomean <- function(x)
{
	x <- x[!is.na(x)]
	ifelse(all(x > 0), exp(mean(log(x))), NA)	# データは全部正の値でなくてはならない。戻り値は，対数値の平均値の逆対数（指数）
}

harmean <- function(x)
{
	x <- x[!is.na(x)]
	ifelse(all(x > 0), 1/mean(1/x), NA)		# データは全部正の値でなくてはならない。戻り値は，逆数の平均値の逆数
}

my.sum <- function(x)					# 大きさの違う数を足し算するときには若干の注意が必要（あまり効果はない）
{
	x <- x[!is.na(x)]
	sum(x[order(abs(x))])
}

skew <- function(x, method = c("Excel", "ordinary"))
{
	method <- match.arg(method)			# 省略可能なパラメータの処理。標準は Excel（SPSS） と同じ計算法
	x <- x[!is.na(x)]				# 欠損値を持つケースを除く
	n <- length(x)					# データ数
	if (method == "Excel") {			# Excel（SPSS）と同じ計算法
		n*my.sum(scale(x)^3)/(n-1)/(n-2)	# scale は元のデータを標準化する関数
	}
	else {
		my.sum(((x-mean(x))/sqrt((n-1)*var(x)/n))^3)/n	# 標準化は分散の平方根を取った標準偏差による
	}
}

kurt <- function(x, method = c("Excel", "ordinary"))
{
	method <- match.arg(method)			# 省略可能なパラメータの処理。標準は Excel（SPSS） と同じ計算法
	x <- x[!is.na(x)]				# 欠損値を持つケースを除く
	n <- length(x)					# データ数
	if (method == "Excel") {			# Excel（SPSS）と同じ計算法
		n*(n+1)*sum(scale(x)^4)/(n-1)/(n-2)/(n-3)-3*(n-1)^2/(n-2)/(n-3)	# scale は元のデータを標準化する関数
	}
	else {
		sum(((x-mean(x))/sqrt((n-1)*var(x)/n))^4)/n-3	# 標準化は分散の平方根を取った標準偏差による
	}
}

stdev <- sd						# 単に名前の違い

stdevp <- function(x)
{
	x <- x[!is.na(x)]				# 欠損値を持つケースを除く
	n <- length(x)					# データ数
	sd(x)*sqrt((n-1)/n)				# 分散の平方根により定義される標準偏差
}

trimmean <- function(x, p)
{
	mean(x, p/2, na.rm=TRUE)
}

varp <- function(x)
{
	x <- x[!is.na(x)]				# 欠損値を持つケースを除く
	n <- length(x)					# データ数
	(n-1)*var(x)/n					# 不偏分散から分散を求める
}

large <- function(x, k)
{
	rev(sort(x[!is.na(x)]))[k]			# 大きい方からk番目の数値
}

small <- function(x, k)
{
	sort(x[!is.na(x)])[k]				# 小さい方からk番目の数値
}

fact <- factorial					# 名前の違い

combin <- choose					# 名前の違い
# ファン・デル・ワーデン検定（permutation.test 関数用）
vdw2.test <- function(	x,			# 第一群のデータベクトル		
			y)			# 第二群のデータベクトル
{
	x <- x[!is.na(x)]			# 欠損値を持つケースを除く
	y <- y[!is.na(y)]			# 欠損値を持つケースを除く
	n1 <- length(x)				# 第一群のサンプルサイズ
	z <- c(x, y)				# データベクトルを一つにまとめる
	n <- length(z)				# 合計したサンプルサイズ
	S <- sum(qnorm((rank(z)[1:n1])/(n+1)))	# 第一群のデータに対する正規化得点の合計
	return(list(statistic=S))
}
# ファン・デル・ワーデン検定
vdw.test <- function(	x,			# 第一群のデータベクトル		
			y)			# 第二群のデータベクトル
{
	method <- "ファン・デル・ワーデン検定"
	data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	x <- x[!is.na(x)]			# 欠損値を持つケースを除く
	y <- y[!is.na(y)]			# 欠損値を持つケースを除く
	n1 <- length(x)				# 第一群のサンプルサイズ
	n2 <- length(y)				# 第二群のサンプルサイズ
	n <- n1+n2				# 合計したサンプルサイズ
	z <- c(x, y)				# データベクトルを一つにまとめる
	z <- qnorm(rank(z)/(n+1))		# 各データに対する正規化得点を計算
	S <- abs(sum(z[1:n1]))			# 第一群のデータに対する正規化得点の合計
	V <- n1*n2/(n^2-n)*sum(z^2)		# 分散
	Z <- S/sqrt(V)				# 検定統計量の正規化得点
	P <- pnorm(Z, lower.tail=FALSE)*2	# P 値
	return(structure(list(statistic=c(S=S, "V(S)"=V, Z=Z), p.value=P,
		method=method, data.name=data.name), class="htest"))
}
version2 <- function() cat("2019/00/01 (17:31:33)")
# 測定値をワイブル確率紙にプロットする
weibull <- function(	x,				# データベクトル
			color="gray")			# 格子線を描く色
{
	weib <- function(p) log10(log10(1/(1-p)))	# ワイブル分布の尺度に変更
	log.axis <- function(z)				# 対数軸を描く関数
	{
		z <- floor(log10(z))			# 対数にしたときの整数部
		log.min <- min(z)			# 最小値
		z2 <- 1:10*10^log.min			# 値の範囲をカバーするように
		n <- max(z)-log.min			# 10 倍しながら順次，右の位置に目盛りを描く
		z2 <- rep(z2, n+1)*10^rep(0:n, each=10)	# 対数目盛り位置の数値
		log.z2 <- log10(z2)			# 目盛りを描く位置
		axis(1, at=log.z2, labels=z2)		# log.z2 の位置に，z2 という数値を描く
		abline(v=log.z2, col=color)		# 垂直格子線を描く
	}

	x <- x[!is.na(x)]				# 欠損値を除く
	n <- length(x)					# 有効データ数
	x <- sort(x)					# 昇順に並べ替える
	log.x <- log10(x)				# 常用対数を取る
	y <- weib(((1:n)-0.5)/n)		# ワイブル分布における累積密度の位置
	y0 <- c(10^(-10:0), 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 99.9, 99.99, 99.999)
	y0 <- y0[y0 > 10/n]
	probs <- weib(y0/100)				# 目盛数値の位置
	plot(c(log.x[1], log.x[n]), c(probs[1], probs[length(probs)]), type="n", xaxt="n", yaxt="n",
		xlab="Observed Value", ylab="Cumulative Percent",
		main="Weibull Probability Paper")
	abline(h=probs, col="grey")			# 水平の格子線
	log.axis(x)					# 横軸を描く
	axis(2, at=probs, labels=y0)			# 縦軸を描く
	points(log.x, y)				#データ点を描く
}
# ワイブル分布のパラメータの最尤推定
weibull.par <- function(x,					# データベクトル
			epsilon=1e-7,				# 収束判定値
			max.loop=500)				# 収束計算の上限回数
{
	x <- x[!is.na(x)]					# 欠損値を持つケースを除く
	n <- length(x)						# データの個数
	m <- a <- 1						# 初期値を設定する
	error <- TRUE						# 収束したかどうかのフラッグ
	for (i in 1:max.loop) {					# 収束するまで繰り返し
		ao <- a						# 改善前の値
		mo <- m						# 改善前の値
		temp1 <- x^mo
		temp2 <- log(x)
		a <- n/sum(temp1)				# 改善後の値
		m <- n/(a*sum(temp1*temp2)-sum(temp2))		# 改善後の値
		if (abs(a-ao) < epsilon &&
		    abs(m-mo) < epsilon) {			# 差が収束判定値以下ならば，
			error <- FALSE				# 収束した
			break
		}
	}
	if (error) {						# 収束しなかった場合 error = TRUE のまま
		warning("収束しませんでした")
	}
	scale <- (1/a)^(1/m)					# 尺度パラメータ
	return(c("shape"=m, "scale"=scale))			# m は形状パラメータ
}
# ウィルコクソンの符号付順位和検定
wilcox.paired.test <- function(	x,		# 対応する測定値ベクトル
				y,		# 対応する測定値ベクトル
				...)		# wilcox.test へ引き渡すその他の引数
{
	ok <- complete.cases(x, y)		# 欠損値を持つケースを除く	
	x <- x[ok]
	y <- y[ok]

	ok = FALSE				# データが全部整数値になったかどうか
	for (i in 1:10) {			# もとのデータを 10 倍ずつ大きくしていって，全部が整数になったら検定を行う
		if (all(c(x, y)*10^i == floor(c(x, y)*10^i))) {
			ok = TRUE
			break
		}
	}
	if (ok == FALSE) {
		stop("Give up!")
	}
	
	x <- x*10^i				# 整数化
	y <- y*10^i
	wilcox.test(x, y, paired=TRUE, ...)	# R が用意している関数を利用する
}
# ウィルコクソン・スコア（平均順位）
wilcoxon.score <- function(n)		# カテゴリーの度数ベクトル
{
	a <- cumsum(n)
	a <- c(0, a[-length(a)])
	return(a+(n+1)/2)		# カテゴリーのウィルコクソン・スコア
}
# 多変量に拡張された平均値の差の検定（ウィルクスのΛ）
wilks <- function(	dat,					# データ行列
			g)					# 群変数
{
	method <- "多変量に拡張された平均値の差の検定（ウィルクスのΛ）"
	data.name <- paste(deparse(substitute(dat)), "~", deparse(substitute(g)))
	OK <- complete.cases(dat, g)				# 欠損値を持つケースを除く
	dat <- dat[OK,]
	g <- as.factor(g[OK])
	nv <- ncol(dat)						# 変数の個数
	g.case <- table(g)					# 各群のケース数
	n <- sum(g.case)					# 全ケース数
	k <- length(g.case)					# 群の個数

	vars <- by(dat, g, function(x) var(x)*(nrow(x)-1))	# 各群の変動・共変動行列
	w <- matrix(rowSums(mapply("+", vars)), nv)		# 群内変動・共変動行列
	t <- (var(dat)*(n-1))					# 全変動・共変動行列

	LAMBDA <- det(w)/det(t)
	sl <- sqrt(LAMBDA)
	if (nv == 2) {						# 2 変数データの場合
		f <- (1-sl)*(n-k-1)/sl/(k-1)
		p <- pf(f, df1 <- 2*(k-1), df2 <- 2*(n-k-1), lower.tail=FALSE)
		result <- list(statistic=c("F value"=f), parameter=c("df1"=df1, "df2"=df2))
	}
	else if (k == 2) {					# 2 群の場合
		f <- (1-LAMBDA)*(n-k-nv+1)/LAMBDA/nv
		p <- pf(f, nv, df2 <- n-k-nv+1, lower.tail=FALSE)
		result <- list(statistic=c("F value"=f), parameter=c("df1"=nv, "df2"=df2))
	}
	else if (k == 3) {					# 3 群の場合
		f <- (1-sl)*(n-k-nv+1)/sl/nv
		p <- pf(f, df1 <- 2*nv, df2 <- 2*(n-k-nv+1), lower.tail=FALSE)
		result <- list(statistic=c("F value"=f), parameter=c("df1"=nv, "df2"=df2))
	}
	else {							# それ以外の場合
		chi.sq <- ((nv+k)/2-n+1)*log(LAMBDA)
		p <- pchisq(chi.sq, df <- nv*(k-1), lower.tail=FALSE)
		result <- list(statistic=c("chi sq. value"=chi.sq), parameter=c("df"=df))
	}
	return(structure(c(result, p.value=p, method=method, data.name=data.name), class="htest"))
}
xtable.prcomp <- function(obj, caption="caption", label="label", pcs=0, digits=3, rev=-1.5, booktabs=FALSE, type=c("latex", "html")) {
# prcomp 関数が返すオブジェクトを LaTeX ソースとして出力する
	loadings <- t(t(obj$rotation)*obj$sdev)
	if (pcs == 0 || pcs > ncol(loadings)) {
		pcs <- sum(colSums(loadings^2) >= 1)
	}
	loadings.output(loadings[, 1:pcs, drop=FALSE], caption, label, digits, rev, booktabs, type=match.arg(type), pc=TRUE, promax=FALSE)
}

xtable.factanal <- function(obj, caption="caption", label="label", digits=3, rev=-1.5, booktabs=FALSE, type=c("latex", "html")) {
# factanal 関数が返すオブジェクトを LaTeX ソースとして出力する
	loadings.output(obj$loadings, caption, label, digits, rev, booktabs, type=match.arg(type), pc=FALSE, promax=any(grepl("promax", obj$call)))
}

loadings.output <- function(loadings, caption, label, digits, rev, booktabs, type, pc, promax) {
	if (pc) {
		pc <- "主成分"
		eig="固有値"
		contr="& 寄与率"
	} else {
		pc <- "因子"
		eig="因子負荷量二乗和"
		contr="& 共通性"
	}
	n <- nrow(loadings)
	factors <- ncol(loadings)
	communality <- rowSums(loadings^2)
	loadings <- cbind(loadings, communality)
	eva <- colSums(loadings^2)
	con <- eva/n*100
	cum <- cumsum(con)
	loadings <- rbind(loadings, eva, con, cum)
	loadings[n+1:3, factors+1] <- NA
	rownames(loadings)[n+1:3] <- c(eig, "寄与率(\\%)", "累積寄与率(\\%)")
	vnames <- rownames(loadings)
#	rotation <- as.character(a$call)
#	print(unlist(strsplit(rotation, " ")))
	if (booktabs) {
		TOPRULE <- "toprule"
		MIDRULE <- "midrule"
		BOTTOMRULE <- "bottomrule"
	} else {
		TOPRULE <- MIDRULE <- BOTTOMRULE <- "hline"
	}
	align <- paste(rep("r", factors+2), collapse="")
	if (type == "latex") {
		cat(sprintf('\\begin{table}[htbp]\n\\caption{%s}\n\\label{%s}\n\\centering\n\\begin{tabular}{%s} \\%s\n', caption, label, align, TOPRULE))
		cat(" &", paste(paste("第", 1:factors, pc, sep=""), collapse=" & "), contr, sprintf("\\\\ \\%s \n", MIDRULE))
		n <- n+3
		for (i in 1:n) {
			cat(vnames[i])
			format <- sprintf(" & $%%.%if$", ifelse(i >= n-1, 1, digits))
			for (j in 1:(factors+1)) {
				if (is.na(loadings[i,j])) {
							cat(" & ")
						}
				else {
					cat(sprintf(format, loadings[i,j]))
				}
			}
			cat("\\\\")
			if (i == n-3) cat(sprintf("\\%s", MIDRULE))
			else if (i < n) cat(sprintf("[%smm]", rev))
			else if (i == n) cat(sprintf("\\%s", BOTTOMRULE))
			cat("\n")
			if (i == n-3 && promax) break
		}
		cat("\\end{tabular}\\end{table}\n")
	}

}

xtable.ftable <- function(	obj,				# crosstabsn が返す ftable オブジェクト
			caption="caption",			# キャプション
			label="label",			# ラベル
			percentage=c("row", "col", "none"),	# % を付ける方向
			same.line=TRUE, 			# % を度数と同じ行に付けるときは TRUE にする
			percentage.font=c("small", "footnotesize", "tiny", "normalsize"), # LaTeX でのフォントサイズの指定 tiny, footnotesize など
			position=c("c", "r", "l"), 		# フィールド内での配置 "c", "r", "l" のいずれか
			rev=-1.5,				# 行間を詰めるための，逆改行の大きさをミリ単位で指定（逆改行しない場合には 0 を指定する）
			booktabs=FALSE)			# TRUE なら \hline の代わりに \toprule, \midrule, \bottomrule を使う
# ftable 関数が返した ftable クラスのオブジェクトを入力し，LaTeX ソースを出力する
# formula で指定する。~ の左辺には1個のみ指定できる
# Sweave から使用するのが便利
# 使用例
# Titanic
# x <- ftable(Survived ~ ., data = Titanic)
# a <- ftable(Survived ~ Sex + Class + Age, data = x)
# xtable(a)
# xtable(ftable(Survived ~ Sex + Class, data = x))

{	row.vars <- attr(obj, "row.vars")
	n.row.vars <- length(row.vars)
	names.row.vars <- names(row.vars)
	m.row.vars <- sapply(row.vars, length)
	col.vars <- attr(obj, "col.vars")
	n.col.vars <- length(col.vars)
	names.col.vars <- names(col.vars)
	m.col.vars <- sapply(col.vars, length)
	if (n.col.vars != 1) {
		stop("col.vars が 1 変数の ftable オブジェクトしか扱えません")
	}
	nrow <- nrow(obj)
	side <- matrix("", nrow, n.row.vars)
	n.block <- nrow/m.row.vars[n.row.vars]
	side[, n.row.vars] <- unlist(rep(row.vars[n.row.vars], n.block))
	for (i in seq_len(n.row.vars-1)) {
		every <- prod(m.row.vars[(i+1):n.row.vars])
		side[(0:(nrow-1))%%every==0, i] <- unlist(row.vars[i])
	}

	percentage <- match.arg(percentage)
	if (percentage == "none") {
		same.line <- FALSE
	}
	percentage.font <- match.arg(percentage.font)
	position <- match.arg(position)

	if (booktabs) {
		toprule <- "\\toprule"
		midrule <- "\\midrule"
	}
	else {
		toprule <- midrule <- "\\hline"
	}

	col.vars <- c(unlist(col.vars[[1]]), "合計")
	fac <- same.line+1
	if (same.line) {
		pos <- c(rep(position, n.row.vars), rep(paste(position, "@{}", position), m.col.vars+1))
		header <- paste(paste(names.row.vars, collapse=" & "), paste("&", paste(col.vars, "\\%", sep=" & ", collapse=" & ")))
		fmt <- sprintf("%%d & {\\%s \\textit{%%6.1f}}", percentage.font)
	}
	else {
		pos <- rep(position, m.col.vars+1+n.row.vars)
		header <- paste(paste(names.row.vars, collapse=" & "), paste(col.vars, collapse=" & "), sep=" & ")
		fmt <- sprintf("{\\%s \\textit{%%5.1f}}", percentage.font)
	}
	cat("\\begin{table}[htbp]\n",
	    "\\caption{", caption, "}\n",
	    "\\label{", label, "}\n",
	    "\\centering\n",
	    "\\begin{tabular}{", pos, "} ", toprule, " \n", sep="")
	cat(paste(rep("&", n.row.vars), collapse=" "))
	cat(sprintf(" \\multicolumn{%i}{c}{%s}\\\\ \\cline{%i-%i}\n",
	    fac*m.col.vars[1], names.col.vars[1], n.row.vars+1, fac*m.col.vars[1]+n.row.vars))
	cat(header, " \\\\ ", midrule, "\n", sep="")

	for (k in 1:n.block) {
		end <- k*m.row.vars[n.row.vars]
		begin <- end-m.row.vars[n.row.vars]+1
		block <- addmargins(obj[begin:end, ])
		side.block <- rbind(side[begin:end, , drop=FALSE ], c(rep("", n.row.vars-1), "合計"))
		if (percentage == "row") {
			pct <- block/block[, m.col.vars+1]*100
		}
		else {
			pct <- t(t(block)/block[m.row.vars[n.row.vars]+1,]*100)
		}
		n <- m.row.vars[n.row.vars]+1
		for (i in 1:n) {
			cat(sprintf("%s &", side.block[i,]))
			if (same.line) {
				cat(gsub("NaN", "---", paste(apply(cbind(block[i,], pct[i,]), 1, function(y) sprintf(fmt, y[1], y[2])), collapse=" & ")))
			}
			else {
				cat(paste(block[i,], collapse=" & "), "\\\\ \n")
				if (percentage != "none") {
					cat(rep(" &", n.row.vars-1))
					cat("\\%", gsub("NaN", "---", sprintf(fmt, pct[i, ])), sep=" & ")
				}
			}
			if (percentage != "none") {
				cat(" \\\\")
			}
			if (i < n-1) {
				cat(sprintf("[%smm]\n", rev))
			}
			else if (i == n) {
				if (end < nrow) {
					cat(sprintf("\\cline{%i-%i}\n", sum(side[end+1,] == "")+1, fac*(m.col.vars[1]+1)+n.row.vars))
				}
				else {
					cat(sprintf("%s\n", toprule))
				}
			}
			else {
				cat(sprintf("\\cline{%i-%i}\n", n.row.vars, fac*(m.col.vars[1]+1)+n.row.vars))
			}
		}
	}
	cat("\\end{tabular}\n",
	    "\\end{table}\n", sep="")
}
xtable.glm <- function(obj, caption="caption", label="label", vif=FALSE, align="lrrrr", digits=rep(3, 4), rev=-1.5, booktabs=FALSE, type=c("latex", "html"), suf=FALSE) {
# glm 関数が返すオブジェクトを LaTeX または html ソースとして出力する
	conv <- function(s) { # 添字を数式モードで
		if (suf) paste0("$", sub("([0-9]+$)", "_{\\1}", s), "$") else s
	}		
	ans2 <- summary(obj)
	ans <- data.frame(ans2$coefficients)
	if (booktabs) {
		TOPRULE <- "toprule"
		MIDRULE <- "midrule"
		BOTTOMRULE <- "bottomrule"
	}
	else {
		TOPRULE <- MIDRULE <- BOTTOMRULE <- "hline"
	}
	ans <- rbind(ans[-1,], ans[1,])
	rownames(ans)[nrow(ans)] <- "定数項"
	colnames(ans) <- c("偏回帰係数", "標準誤差", "$z$値", "$P$値")
	n <- nrow(ans)
	if (match.arg(type) == "latex") {
		cat(sprintf('\\begin{table}[htbp]\n\\caption{%s}\n\\label{%s}\n\\centering\n\\begin{tabular}{%s} \\%s\n', caption, label, align, TOPRULE))
	
		cat(paste(c("", colnames(ans)), collapse=" & "))
		cat(sprintf("\\\\ \\%s \n", MIDRULE))
		for (i in 1:n) {
			cat(conv(rownames(ans)[i]))
			for (j in 1:4) {
				if (is.na(ans[i,j])) {
					cat(" & ")
				}
				else if (j == 4 && ans[i,j] < 0.001) {
					cat(" & $< 0.001$")
				}
				else {
					format <- sprintf(" & $%%.%if$", digits[j])
					cat(sprintf(format, ans[i,j]))
				}
			}
			cat("\\\\")
			if (i < n-1) cat(sprintf("[%smm]", rev))
			if (i == n-1) cat(sprintf("\\%s\n", MIDRULE))
			cat("\n")
		}
		cat(sprintf("\\%s\n", BOTTOMRULE))
		cat(" \\end{tabular}\\end{table}\n")
	}
	else {
		align <- unlist(strsplit(align, ""))[-1]
		align <- sub("r", "right", align)
		align <- sub("l", "left", align)
		align <- sub("c", "center", align)
		colnames(ans) <- c("偏回帰係数", "標準誤差", "\\(z\\) 値", "\\(P\\) 値")
		cat("<TABLE border=1>\n")
		cat(sprintf("<CAPTION ALIGN='top'> %s </CAPTION>\n", caption))
		cat(paste(c("<TR> <TH> ", colnames(ans)), collapse=" </TH> <TH> "))
		cat(" </TH> </TR>\n")
		for (i in 1:n) {
			cat("  <TR> <TD>", rownames(ans)[i])
			for (j in 1:4) {
				if (is.na(ans[i,j])) {
					cat(sprintf(" </TD> <TD align='%s'> ", align[j]))
				}
				else if (j == 4 && ans[i,j] < 0.001) {
					cat(sprintf(" </TD> <TD align='%s'> &lt; 0.001 ", align[j]))
				}
				else {
					format <- sprintf(" </TD> <TD align='%s'> %%.%if", align[j], digits[j])
					cat(sprintf(format, ans[i,j]))
				}
			}
			cat(" </TD> </TR>\n")
		 }
		cat("	</TABLE>\n")
	}
}
xtable.lm <- function(obj, caption="caption", label="label", vif=FALSE, align="lrrrrrr", digits=rep(3, 6), footnote=rep(TRUE, 4), rev=-1.5, booktabs=FALSE, type=c("latex", "html"), suf=FALSE) {
# lm 関数が返すオブジェクトを LaTeX または html ソースとして出力する
	conv <- function(s) { # 添字を数式モードで
		if (suf) paste0("$", sub("([0-9]+$)", "_{\\1}", s), "$") else s
	}		
	ans2 <- summary(obj)
	ans <- data.frame(ans2$coefficients)
	df <- obj$model
	d <- model.matrix(obj$terms, eval(obj$model, parent.frame()))
	ans$std.est <- c(NA, obj$coefficients[-1] * apply(d[, -1, drop=FALSE], 2, sd) / sd(obj$model[, 1]))
	ans$append <- c(NA, diag(solve(cor(d[,-1, drop=FALSE]))))
	name.append <- "VIF"
	if (!vif) {
		ans$append <- 1/ans$append
		name.append <- "トレランス"
	}
	if (booktabs) {
		TOPRULE <- "toprule"
		MIDRULE <- "midrule"
		BOTTOMRULE <- "bottomrule"
	}
	else {
		TOPRULE <- MIDRULE <- BOTTOMRULE <- "hline"
	}
	ans <- rbind(ans[-1,], ans[1,])
	rownames(ans)[nrow(ans)] <- "定数項"
	colnames(ans) <- c("偏回帰係数", "標準誤差", "$t$値", "$P$値", "標準化偏回帰係数", name.append)
	n <- nrow(ans)
	f <- ans2$fstatistic[1]
	df1 <- ans2$fstatistic[2]
	df2 <- ans2$fstatistic[3]
	p <-  pf(f, df1, df2, lower.tail=FALSE)
	if (p < 0.001) {
		p <- " < 0.001"
	} else {
		p <- sprintf(" = %.3f", p)
	}
	if (match.arg(type) == "latex") {
		cat(sprintf('\\begin{table}[htbp]\n\\caption{%s}\n\\label{%s}\n\\centering\n\\begin{tabular}{%s} \\%s\n', caption, label, align, TOPRULE))
	
		cat(paste(c("", colnames(ans)), collapse=" & "))
		cat(sprintf("\\\\ \\%s \n", MIDRULE))
		for (i in 1:n) {
			cat(conv(rownames(ans)[i]))
			for (j in 1:6) {
				if (is.na(ans[i,j])) {
					cat(" & ")
				}
				else if (j == 4 && ans[i,j] < 0.001) {
					cat(" & $< 0.001$")
				}
				else {
					format <- sprintf(" & $%%.%if$", digits[j])
					cat(sprintf(format, ans[i,j]))
				}
			}
			cat("\\\\")
			if (i < n-1) cat(sprintf("[%smm]", rev))
			if (i == n-1) cat(sprintf("\\%s\n", MIDRULE))
			cat("\n")
		}
		cat(sprintf("\\%s\n", BOTTOMRULE))
		if (footnote[1]) {
			cat(sprintf("&\\multicolumn{6}{l}{重相関係数 $R = %.3f$}\\\\[%smm]\n", sqrt(ans2$r.squared), rev))
		}
		if (footnote[2]) {
			cat(sprintf("&\\multicolumn{6}{l}{重相関係数の二乗（決定係数）$R^2 = %.3f$}\\\\[%smm]\n", ans2$r.squared, rev))
		}
		if (footnote[3]) {
			cat(sprintf("&\\multicolumn{6}{l}{自由度調整済み重相関係数の二乗 $= %.3f$}\\\\[%smm]\n", ans2$adj.r.squared, rev))
		}
		if (footnote[4]) {
			cat(sprintf("&\\multicolumn{6}{l}{回帰の分散分析：$F値(%i, %i) = %.3f$，$P 値%s$}\\\\[%smm]\n", df1, df2, f, p, rev))
		}
		cat(" \\end{tabular}\\end{table}\n")
	}
	else {
		align <- unlist(strsplit(align, ""))[-1]
		align <- sub("r", "right", align)
		align <- sub("l", "left", align)
		align <- sub("c", "center", align)
		colnames(ans) <- c("偏回帰係数", "標準誤差", "\\(t\\) 値", "\\(P\\) 値", "標準化偏回帰係数", name.append)
		cat("<TABLE border=1>\n")
		cat(sprintf("<CAPTION ALIGN='top'> %s </CAPTION>\n", caption))
		cat(paste(c("<TR> <TH> ", colnames(ans)), collapse=" </TH> <TH> "))
		cat(" </TH> </TR>\n")
		for (i in 1:n) {
			cat("  <TR> <TD>", rownames(ans)[i])
			for (j in 1:6) {
				if (is.na(ans[i,j])) {
					cat(sprintf(" </TD> <TD align='%s'> ", align[j]))
				}
				else if (j == 4 && ans[i,j] < 0.001) {
					cat(sprintf(" </TD> <TD align='%s'> &lt; 0.001 ", align[j]))
				}
				else {
					format <- sprintf(" </TD> <TD align='%s'> %%.%if", align[j], digits[j])
					cat(sprintf(format, ans[i,j]))
				}
			}
			cat(" </TD> </TR>\n")
		 }
		if (footnote[1]) {
			cat(sprintf("<TR> <TD colspan=7> 重相関係数 \\(R\\) = %.3f </TD> </TR>\n", sqrt(ans2$r.squared), rev))
		}
		if (footnote[2]) {
			cat(sprintf("<TR> <TD colspan=7> 重相関係数の二乗（決定係数）\\(R^2\\) = %.3f </TD> </TR>\n", ans2$r.squared, rev))
		}
		if (footnote[3]) {
			cat(sprintf("<TR> <TD colspan=7> 自由度調整済み重相関係数の二乗 = %.3f </TD> </TR>\n", ans2$adj.r.squared, rev))
		}
		if (footnote[4]) {
        	if (p == " < 0.001") {
        		p <- " &lt; 0.001"
        	}
			cat(sprintf("<TR> <TD colspan=7> 回帰の分散分析：\\(F\\)値(%i, %i) = %.3f，\\(P\\) 値%s </TD> </TR>\n", df1, df2, f, p))
		}
		cat("	</TABLE>\n")
	}
}