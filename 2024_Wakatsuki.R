library(openxlsx)



#### 1. パスの設定 ####

# 結果の保存先
path_out <- file.path(
  '..', '2-2_Paper', '2024_Wakatsuki', 'CleanData', 
  paste0(Sys.Date(), '_Mutsuhomare'))
dir.create(path_out, recursive = TRUE)

# 読み込むデータ
path_dat <- file.path(
  '..', '2-2_Paper', '2024_Wakatsuki', 'RawData', 'CG_database.xlsx')



#### 2. データの読み込み ####

rawdat <- list()

# 列名を読み込み
rawdat$head <- read.xlsx(path_dat, rows = 1, colNames = FALSE)

# データ本体を読み込んで列名を付与
rawdat$all <- read.xlsx(path_dat, startRow = 3, colNames = FALSE, na.strings = 'na')
colnames(rawdat$all) <- unlist(rawdat$head)
# View(rawdat$all)



#### 3. データの形式の整形 ####

# 列名を置き換え
dat <- list(df = rawdat$all)
colnames(dat$df)[colnames(dat$df) == 'Pref.']                <- 'Pref'
colnames(dat$df)[colnames(dat$df) == 'Institution/Station']  <- 'InstOrStation'
colnames(dat$df)[colnames(dat$df) == 'Soil Groups ']         <- 'SoilGroups'
colnames(dat$df)[colnames(dat$df) == 'Cultivar name (E)']    <- 'Cultivar'
colnames(dat$df)[colnames(dat$df) == 'DOY']                  <- 'HeadingDOY'
colnames(dat$df)[colnames(dat$df) == 'Head rice']            <- 'HeadRice'
colnames(dat$df)[colnames(dat$df) == 'Milky-white']          <- 'MilkyWhite'
colnames(dat$df)[colnames(dat$df) == 'Basal-white']          <- 'BasalWhite'
colnames(dat$df)[colnames(dat$df) == 'Belly white']          <- 'BellyWhite'
colnames(dat$df)[colnames(dat$df) == 'sieve']                <- 'Sieve'
colnames(dat$df)[colnames(dat$df) == 'Measurement methods '] <- 'MeasurementMethods'
colnames(dat$df)[colnames(dat$df) == 'URL/doi']              <- 'URLorDOI'
colnames(dat$df)[colnames(dat$df) == 'Ref_No.']              <- 'RefNo'

# データフレーム中の値に、末尾にスペースが入っているものがあったら、
# それを取り除く
remove_tail_space <- function(dat) {
  
  for (cn in colnames(dat)) {
    val <- dat[, cn]
    val_with_space <- unique(val[grepl(' $', val)])
    if (length(val_with_space) >= 1) {
      cat(sprintf('Removing unnecessary space in "%s".\n', cn))
      for (txt in val_with_space) {
        cat(sprintf('  "%s"\n', txt))
      }
      dat[, cn] <- gsub(' $', '', val)
    }
  }
  dat
  
}
dat$df <- remove_tail_space(dat$df)



#### 4. データの中身の整形 ####

# 列の情報を抽出
dat$col_info <- data.frame(name = colnames(dat$df))
dat$col_info <- within(dat$col_info, {
  
  # データの中でどの列に値が入っているか登録
  is_val <- seq_along(name) %in% 12:25
  
  # 値の種類を登録
  val_type <- ifelse(
    grepl('^(Tave)|(TaHD..)|(SR)|(RH)$', name), 'Weather', 'Phenotype')
  val_type[!is_val] <- NA
  
})

# 異常値を検出
detect_strange_value <- function(df, colinfo) {
  
  # 単位が%のデータを抽出
  is_pheno <- rep(FALSE, nrow(colinfo))
  is_pheno[colinfo$val_type == 'Phenotype'] <- TRUE
  is_pheno[colinfo$name == 'HeadingDOY'] <- FALSE
  pheno_df <- df[, is_pheno]
  
  # ％の表現型が1~99に収まっていないものを見つける
  is_strange <- pheno_df >= 100 | pheno_df <= 0
  is_strange[is.na(is_strange)] <- FALSE
  
  # 見つかったデータを表示
  cat('value\n')
  print(pheno_df[rowSums(is_strange) >= 1, ])
  cat('\nis_strange\n')
  print(is_strange[rowSums(is_strange) >= 1, ])
  
  pheno_df[is_strange] <- NA
  df[, is_pheno] <- pheno_df
  df
  
}
dat_clean <- dat
dat_clean$df <- detect_strange_value(dat_clean$df, dat_clean$col_info)

# 品種名を修正する
dat_clean$df$Cultivar <- gsub('Mutsuhomarre', 'Mutsuhomare', dat_clean$df$Cultivar)



#### 5. 情報の整理 ####

# データを情報と値に分解
dat_clean <- within(dat_clean, {
  mat <- as.matrix(df[, col_info$is_val])
  row_info <- df[, !col_info$is_val]
})

# ReferenceとSiteの情報を抽出
ref_site_info <- dat_clean$row_info[, c('Pref', 'Reference', 'RefNo')]
ref_site_info <- unique(ref_site_info)
ref_site_info <- ref_site_info[order(ref_site_info$RefNo), ]

# 解析に必要なreference-site IDを作成
create_ref_site_id <- function(ref_site_info) {
  
  # Referenceのリストを作成
  ref_info <- unique(ref_site_info[, c('Reference', 'RefNo')])
  
  # Reference名の略称を作成
  ref_info$RefAbb <- ref_info$Reference
  ref_info$RefAbb <- gsub('( et? al\\.?)?, ', '', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('^RCC database- "', 'RCC-', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('"-?$', '', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('The Performance Tests for Recommended Varieties of Rice database', 
                          'RcommendedVarietyTest', ref_info$RefAbb)
  
  # Referenceの略称が重複しているものについて、末尾に小文字のアルファベットを付与
  ref_abb_dup <- unique(ref_info$RefAbb[duplicated(ref_info$RefAbb)])
  for (refabb0 in ref_abb_dup) {
    # refabb0 <- ref_abb_dup[1]
    r <- which(ref_info$RefAbb == refabb0)
    ref_info$RefAbb[r] <- paste0(ref_info$RefAbb[r], letters[seq_along(r)])
  }
  
  # Referenceの略称を、ref_site_infoにくっつける
  ref_site_info <- merge(ref_site_info, ref_info, all = TRUE, sort = FALSE)
  
  # Referenceの略称とsiteの情報を繋げ、IDとする
  ref_site_info$RefSiteID <- with(ref_site_info, paste0(RefAbb, '_', Pref))
  ref_site_info$RefAbb <- NULL
  return(ref_site_info)
  
}
ref_site_info <- create_ref_site_id(ref_site_info)

# 作成したIDをデータに付与
add_column_of_ref_site_id <- function(ref_site_info, info) {
  # データフレームinfoに、Reference-site IDを付与する
  # ref_site_infoは、ReferenceとSiteの組の一覧表
  
  # 万が一、一覧表のIDに重複があったら、エラーとして処理を停止
  if (any(duplicated(ref_site_info$RefSiteID))) {
    stop('Given RefSiteID is duplicated.')
  }
  
  # Site-Refの一覧表を、IDとそれ以外に切り離す
  ref_site_id <- ref_site_info$RefSiteID
  ref_site_info <- ref_site_info[, colnames(ref_site_info) != 'RefSiteID']
  
  # infoから必要な情報だけ取り出す
  info0 <- info[, colnames(ref_site_info)]
  
  # IDを付与
  id <- vector('numeric', nrow(info))
  colsole_width <- getOption('width')
  cat('Searching reference-site ID...\n')
  for (i in seq_along(id)) {
    
    info0_now <- info0[i, ]
    for (j in seq_along(ref_site_id)) {
      
      # 一覧表の中に一致する情報が見つかったら、IDを付与するためループを抜ける
      if (all(info0_now == ref_site_info[j, ])) break
      
      # 万が一、一覧表の中に対応する情報が見つからなかったら、
      # エラーとして処理を停止
      if (j == length(ref_site_id)) {
        cat('No matching information was found for:\n')
        cat(paste(rep('-', colsole_width), collapse = ''))
        cat('\n')
        print(info0_now)
        cat(paste(rep('-', colsole_width), collapse = ''))
        cat('\n')
        stop()
      }
      
    }
    id[i] <- ref_site_id[j]
    
  }
  info$RefSiteID <- id
  cat('Finished.\n')
  return(info)
  
}
dat_clean$row_info <- add_column_of_ref_site_id(ref_site_info, dat_clean$row_info)

# ReferenceとSiteの組合せに番号を振る
ref_site_info$RefSiteNo <- seq_len(nrow(ref_site_info))
dat_clean$row_info$RefSiteNo <- with(
  ref_site_info, RefSiteNo[match(dat_clean$row_info$RefSiteID, RefSiteID)])

# 結果を保存
saveRDS(dat_clean, file.path(path_out, 'Data_Main.rds'))
saveRDS(ref_site_info, file.path(path_out, 'Info_Reference-Site.rds'))




#### 6. 気象要因をセンタリング ####

dat_center <- within(dat_clean, {
  
  # 気象要因の平均を計算
  is_weather <- rep(FALSE, nrow(col_info))
  is_weather[col_info$val_type == 'Weather'] <- TRUE
  weather_mean <- colMeans(dat_clean$df[, is_weather])
  
  # 気象要因の平均をcol_infoの中に記録
  col_info$mean <- NA
  col_info$mean[is_weather] <- weather_mean
  
  # 気象要因をセンタリング
  df[, is_weather] <- df[, is_weather] - rep(weather_mean, each = nrow(df))
  mat <- df[, col_info$is_val]
  rm(list = c('is_weather', 'weather_mean'))
  
})
saveRDS(dat_center, file.path(path_out, 'Data_Center.rds'))



#### 7. 値をチェック ####

for (data_name in c('Main', 'Center')) {
  
  dat_now <- if (data_name == 'Main') dat_clean else dat_center
  pdf(file.path(path_out, sprintf('Fig_%s.pdf', data_name)),
      height = 7, width = 9)
  for (i in seq_len(ncol(dat_now$mat))) {
    par(mar = c(8, 4, 2, 2))
    b <- with(dat_now, boxplot(
      mat[, i] ~ row_info$RefSiteID, main = colnames(mat)[i], 
      xlab = '', ylab = colnames(mat)[i], xaxt = 'n'))
    abline(v = seq_along(b$names), col = rgb(0, 0, 0, 0.2))
    text(x = seq_along(b$names), 
         y = sum(par('usr')[3:4] * c(1.01, -0.01)),
         labels = b$names, cex = 0.6,
         adj = c(1, 0.5), srt = 90, xpd = TRUE)
  }
  dev.off()
  
}
