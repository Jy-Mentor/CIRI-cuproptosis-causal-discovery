# === 通过 multiMiR R 包提取 miRTarBase 数据 ===
library(multiMiR)

output_dir <- "D:/反向网络药理学/GAT拓展维度"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("multiMiR 版本:", as.character(packageVersion("multiMiR")), "\n")

# 测试小批次查询
test_mirnas <- c("hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p")
cat("测试查询 (", paste(test_mirnas, collapse = ", "), ")\n", sep = "")

result <- get_multimir(
  mirna = test_mirnas,
  org = "hsa",
  table = "validated",
  predicted.cutoff.type = "n",
  use.tibble = FALSE,
  legacy.out = FALSE
)

cat("测试结果行数:", nrow(result@data), "\n")
cat("列名:", paste(colnames(result@data), collapse = ", "), "\n")

# 批量查询所有已知人类miRNA
all_mirnas <- unique(result@data[[grep("mirna", colnames(result@data), ignore.case=TRUE, value=TRUE)[1]]])
cat("从测试结果获取的miRNA数:", length(all_mirnas), "\n")

# 手动补充已知的人类miRNA主要列表
major_mirnas <- c(
  "hsa-let-7a-5p", "hsa-let-7b-5p", "hsa-let-7c-5p", "hsa-let-7d-5p", "hsa-let-7e-5p",
  "hsa-let-7f-5p", "hsa-let-7g-5p", "hsa-let-7i-5p",
  "hsa-miR-1-3p", "hsa-miR-7-5p", "hsa-miR-9-5p", "hsa-miR-10a-5p", "hsa-miR-10b-5p",
  "hsa-miR-15a-5p", "hsa-miR-15b-5p", "hsa-miR-16-5p", "hsa-miR-17-5p", "hsa-miR-18a-5p",
  "hsa-miR-19a-3p", "hsa-miR-19b-3p", "hsa-miR-20a-5p", "hsa-miR-20b-5p",
  "hsa-miR-21-5p", "hsa-miR-22-3p", "hsa-miR-23a-3p", "hsa-miR-23b-3p", "hsa-miR-24-3p",
  "hsa-miR-25-3p", "hsa-miR-26a-5p", "hsa-miR-26b-5p", "hsa-miR-27a-3p", "hsa-miR-27b-3p",
  "hsa-miR-28-5p", "hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p",
  "hsa-miR-30a-5p", "hsa-miR-30b-5p", "hsa-miR-30c-5p", "hsa-miR-30d-5p", "hsa-miR-30e-5p",
  "hsa-miR-31-5p", "hsa-miR-32-5p", "hsa-miR-33a-5p", "hsa-miR-33b-5p",
  "hsa-miR-34a-5p", "hsa-miR-34b-3p", "hsa-miR-34c-5p",
  "hsa-miR-92a-3p", "hsa-miR-93-5p", "hsa-miR-96-5p",
  "hsa-miR-100-5p", "hsa-miR-101-3p", "hsa-miR-103a-3p", "hsa-miR-105-5p",
  "hsa-miR-106a-5p", "hsa-miR-106b-5p", "hsa-miR-107",
  "hsa-miR-122-5p", "hsa-miR-124-3p", "hsa-miR-125a-5p", "hsa-miR-125b-5p",
  "hsa-miR-126-3p", "hsa-miR-127-3p", "hsa-miR-128-3p",
  "hsa-miR-130a-3p", "hsa-miR-130b-3p", "hsa-miR-132-3p", "hsa-miR-133a-3p", "hsa-miR-133b",
  "hsa-miR-134-5p", "hsa-miR-135a-5p", "hsa-miR-135b-5p",
  "hsa-miR-136-5p", "hsa-miR-137", "hsa-miR-138-5p", "hsa-miR-139-5p",
  "hsa-miR-140-3p", "hsa-miR-141-3p", "hsa-miR-142-3p", "hsa-miR-143-3p", "hsa-miR-144-3p",
  "hsa-miR-145-5p", "hsa-miR-146a-5p", "hsa-miR-146b-5p",
  "hsa-miR-147a", "hsa-miR-148a-3p", "hsa-miR-148b-3p",
  "hsa-miR-149-5p", "hsa-miR-150-5p", "hsa-miR-151a-3p",
  "hsa-miR-152-3p", "hsa-miR-153-3p", "hsa-miR-154-5p",
  "hsa-miR-155-5p", "hsa-miR-181a-5p", "hsa-miR-181b-5p", "hsa-miR-181c-5p", "hsa-miR-181d-5p",
  "hsa-miR-182-5p", "hsa-miR-183-5p", "hsa-miR-184",
  "hsa-miR-185-5p", "hsa-miR-186-5p", "hsa-miR-187-3p",
  "hsa-miR-188-5p", "hsa-miR-190a-5p", "hsa-miR-191-5p",
  "hsa-miR-192-5p", "hsa-miR-193a-3p", "hsa-miR-193b-3p",
  "hsa-miR-194-5p", "hsa-miR-195-5p", "hsa-miR-196a-5p", "hsa-miR-196b-5p",
  "hsa-miR-197-3p", "hsa-miR-198", "hsa-miR-199a-3p", "hsa-miR-199a-5p", "hsa-miR-199b-5p",
  "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p",
  "hsa-miR-203a-3p", "hsa-miR-204-5p", "hsa-miR-205-5p",
  "hsa-miR-206", "hsa-miR-208a-3p", "hsa-miR-208b-3p",
  "hsa-miR-210-3p", "hsa-miR-211-5p", "hsa-miR-212-3p",
  "hsa-miR-214-3p", "hsa-miR-215-5p", "hsa-miR-216a-5p",
  "hsa-miR-217", "hsa-miR-218-5p", "hsa-miR-219a-5p",
  "hsa-miR-221-3p", "hsa-miR-222-3p", "hsa-miR-223-3p",
  "hsa-miR-224-5p", "hsa-miR-296-5p", "hsa-miR-299-3p",
  "hsa-miR-301a-3p", "hsa-miR-301b-3p", "hsa-miR-302a-3p",
  "hsa-miR-320a", "hsa-miR-323a-3p", "hsa-miR-324-5p",
  "hsa-miR-326", "hsa-miR-328-3p", "hsa-miR-329-3p",
  "hsa-miR-330-5p", "hsa-miR-331-3p", "hsa-miR-335-5p",
  "hsa-miR-337-3p", "hsa-miR-338-3p", "hsa-miR-339-5p",
  "hsa-miR-340-5p", "hsa-miR-342-3p", "hsa-miR-345-5p",
  "hsa-miR-346", "hsa-miR-361-5p", "hsa-miR-362-5p",
  "hsa-miR-363-3p", "hsa-miR-365a-3p", "hsa-miR-365b-3p",
  "hsa-miR-370-3p", "hsa-miR-371a-5p", "hsa-miR-372-3p",
  "hsa-miR-373-3p", "hsa-miR-374a-5p", "hsa-miR-374b-5p",
  "hsa-miR-375", "hsa-miR-376a-3p", "hsa-miR-376c-3p",
  "hsa-miR-377-3p", "hsa-miR-378a-3p", "hsa-miR-379-5p",
  "hsa-miR-380-5p", "hsa-miR-381-3p", "hsa-miR-382-5p",
  "hsa-miR-383-5p", "hsa-miR-384", "hsa-miR-409-3p",
  "hsa-miR-410-3p", "hsa-miR-411-5p", "hsa-miR-412-5p",
  "hsa-miR-421", "hsa-miR-422a", "hsa-miR-423-5p",
  "hsa-miR-424-5p", "hsa-miR-425-5p", "hsa-miR-429",
  "hsa-miR-431-5p", "hsa-miR-432-5p", "hsa-miR-433-3p",
  "hsa-miR-449a", "hsa-miR-449b-5p", "hsa-miR-450a-5p",
  "hsa-miR-451a", "hsa-miR-452-5p", "hsa-miR-454-3p",
  "hsa-miR-455-5p", "hsa-miR-483-5p", "hsa-miR-484",
  "hsa-miR-485-5p", "hsa-miR-486-5p", "hsa-miR-487a-3p",
  "hsa-miR-488-3p", "hsa-miR-489-3p", "hsa-miR-490-3p",
  "hsa-miR-491-5p", "hsa-miR-493-5p", "hsa-miR-494-3p",
  "hsa-miR-495-3p", "hsa-miR-496", "hsa-miR-497-5p",
  "hsa-miR-498", "hsa-miR-499a-5p", "hsa-miR-500a-5p",
  "hsa-miR-501-5p", "hsa-miR-502-5p", "hsa-miR-503-5p",
  "hsa-miR-504-5p", "hsa-miR-505-5p", "hsa-miR-506-3p",
  "hsa-miR-507", "hsa-miR-508-3p", "hsa-miR-509-3p",
  "hsa-miR-510", "hsa-miR-511", "hsa-miR-512-5p",
  "hsa-miR-513a-5p", "hsa-miR-514a-3p", "hsa-miR-515-5p",
  "hsa-miR-516a-5p", "hsa-miR-517a-3p", "hsa-miR-518a-5p",
  "hsa-miR-518b", "hsa-miR-518c-5p", "hsa-miR-518d-5p",
  "hsa-miR-518e-5p", "hsa-miR-518f-5p", "hsa-miR-519a-5p",
  "hsa-miR-519b-5p", "hsa-miR-519c-3p", "hsa-miR-519d",
  "hsa-miR-520a-3p", "hsa-miR-520b", "hsa-miR-520c-3p",
  "hsa-miR-520d-3p", "hsa-miR-520e", "hsa-miR-520f-3p",
  "hsa-miR-520g", "hsa-miR-520h", "hsa-miR-521",
  "hsa-miR-522-3p", "hsa-miR-523-3p", "hsa-miR-524-5p",
  "hsa-miR-525-5p", "hsa-miR-526b-5p", "hsa-miR-527",
  "hsa-miR-532-5p", "hsa-miR-539-5p", "hsa-miR-541-3p",
  "hsa-miR-542-5p", "hsa-miR-543", "hsa-miR-544a",
  "hsa-miR-545-5p", "hsa-miR-548a-3p", "hsa-miR-548b-3p",
  "hsa-miR-548c-3p", "hsa-miR-548d-3p", "hsa-miR-548e-3p",
  "hsa-miR-548f-3p", "hsa-miR-548g-3p", "hsa-miR-548h-5p",
  "hsa-miR-548i", "hsa-miR-548j-5p", "hsa-miR-548k",
  "hsa-miR-548l", "hsa-miR-548m", "hsa-miR-548n",
  "hsa-miR-548o", "hsa-miR-548p", "hsa-miR-548q",
  "hsa-miR-548s", "hsa-miR-548t-5p", "hsa-miR-548u",
  "hsa-miR-548v", "hsa-miR-548w", "hsa-miR-548x-3p",
  "hsa-miR-548y", "hsa-miR-549a", "hsa-miR-550a-3p",
  "hsa-miR-551a", "hsa-miR-551b-3p", "hsa-miR-552-3p",
  "hsa-miR-553", "hsa-miR-554", "hsa-miR-555",
  "hsa-miR-556-3p", "hsa-miR-557", "hsa-miR-558",
  "hsa-miR-559", "hsa-miR-560", "hsa-miR-561-3p",
  "hsa-miR-562", "hsa-miR-563", "hsa-miR-564",
  "hsa-miR-565", "hsa-miR-566", "hsa-miR-567",
  "hsa-miR-568", "hsa-miR-569", "hsa-miR-570-3p",
  "hsa-miR-571", "hsa-miR-572", "hsa-miR-573",
  "hsa-miR-574-3p", "hsa-miR-575", "hsa-miR-576-5p",
  "hsa-miR-577", "hsa-miR-578", "hsa-miR-579-3p",
  "hsa-miR-580-3p", "hsa-miR-581", "hsa-miR-582-3p",
  "hsa-miR-583", "hsa-miR-584-5p", "hsa-miR-585-3p",
  "hsa-miR-586", "hsa-miR-587", "hsa-miR-588",
  "hsa-miR-589-3p", "hsa-miR-590-3p", "hsa-miR-591",
  "hsa-miR-592", "hsa-miR-593-5p", "hsa-miR-595",
  "hsa-miR-596", "hsa-miR-597", "hsa-miR-598-3p",
  "hsa-miR-599", "hsa-miR-600", "hsa-miR-601",
  "hsa-miR-602", "hsa-miR-603", "hsa-miR-604",
  "hsa-miR-605-5p", "hsa-miR-606", "hsa-miR-606",
  "hsa-miR-607", "hsa-miR-608", "hsa-miR-609",
  "hsa-miR-610", "hsa-miR-611", "hsa-miR-612",
  "hsa-miR-613", "hsa-miR-614", "hsa-miR-615-5p",
  "hsa-miR-616-3p", "hsa-miR-617", "hsa-miR-618",
  "hsa-miR-619-5p", "hsa-miR-620", "hsa-miR-621",
  "hsa-miR-622", "hsa-miR-623", "hsa-miR-624-5p",
  "hsa-miR-625-5p", "hsa-miR-626", "hsa-miR-627-5p",
  "hsa-miR-628-5p", "hsa-miR-629-5p", "hsa-miR-630",
  "hsa-miR-631", "hsa-miR-632", "hsa-miR-633",
  "hsa-miR-634", "hsa-miR-635", "hsa-miR-636",
  "hsa-miR-637", "hsa-miR-638", "hsa-miR-639",
  "hsa-miR-640", "hsa-miR-641", "hsa-miR-642a-5p",
  "hsa-miR-643", "hsa-miR-644a", "hsa-miR-645",
  "hsa-miR-646", "hsa-miR-647", "hsa-miR-648",
  "hsa-miR-649", "hsa-miR-650", "hsa-miR-651-5p",
  "hsa-miR-652-3p", "hsa-miR-653-5p", "hsa-miR-654-3p",
  "hsa-miR-655-3p", "hsa-miR-656-3p", "hsa-miR-657",
  "hsa-miR-658", "hsa-miR-659-3p", "hsa-miR-660-5p",
  "hsa-miR-661", "hsa-miR-662", "hsa-miR-663a",
  "hsa-miR-664a-3p", "hsa-miR-665", "hsa-miR-668-3p",
  "hsa-miR-670-3p", "hsa-miR-671-5p", "hsa-miR-673-5p",
  "hsa-miR-674-5p", "hsa-miR-675-5p", "hsa-miR-676-3p",
  "hsa-miR-758-3p", "hsa-miR-759", "hsa-miR-760",
  "hsa-miR-761", "hsa-miR-762", "hsa-miR-764",
  "hsa-miR-765", "hsa-miR-766-3p", "hsa-miR-767-5p",
  "hsa-miR-768-3p", "hsa-miR-769-5p", "hsa-miR-770-5p",
  "hsa-miR-802", "hsa-miR-873-5p", "hsa-miR-874-3p",
  "hsa-miR-875-5p", "hsa-miR-876-5p", "hsa-miR-877-5p",
  "hsa-miR-885-5p", "hsa-miR-886-3p", "hsa-miR-887-3p",
  "hsa-miR-888-5p", "hsa-miR-889-3p", "hsa-miR-890",
  "hsa-miR-891a-5p", "hsa-miR-891b", "hsa-miR-892a",
  "hsa-miR-892b", "hsa-miR-920", "hsa-miR-921",
  "hsa-miR-922", "hsa-miR-923", "hsa-miR-924",
  "hsa-miR-925", "hsa-miR-933", "hsa-miR-934",
  "hsa-miR-935", "hsa-miR-936", "hsa-miR-937-5p",
  "hsa-miR-938", "hsa-miR-939-5p", "hsa-miR-940",
  "hsa-miR-941", "hsa-miR-942-5p", "hsa-miR-943",
  "hsa-miR-944", "hsa-miR-1200", "hsa-miR-1201",
  "hsa-miR-1202", "hsa-miR-1203", "hsa-miR-1204",
  "hsa-miR-1205", "hsa-miR-1206", "hsa-miR-1207-5p",
  "hsa-miR-1208", "hsa-miR-1224-5p", "hsa-miR-1225-5p",
  "hsa-miR-1226-3p", "hsa-miR-1227-3p", "hsa-miR-1228-3p",
  "hsa-miR-1229-3p", "hsa-miR-1231", "hsa-miR-1233-3p",
  "hsa-miR-1234-3p", "hsa-miR-1236-3p", "hsa-miR-1237-3p",
  "hsa-miR-1238-3p", "hsa-miR-1243", "hsa-miR-1244",
  "hsa-miR-1245a", "hsa-miR-1246", "hsa-miR-1247-5p",
  "hsa-miR-1248", "hsa-miR-1249-5p", "hsa-miR-1250-5p",
  "hsa-miR-1251-5p", "hsa-miR-1252-5p", "hsa-miR-1253",
  "hsa-miR-1254", "hsa-miR-1255a", "hsa-miR-1255b-5p",
  "hsa-miR-1256", "hsa-miR-1257", "hsa-miR-1258",
  "hsa-miR-1259", "hsa-miR-1260a", "hsa-miR-1260b",
  "hsa-miR-1261", "hsa-miR-1262", "hsa-miR-1263",
  "hsa-miR-1264", "hsa-miR-1265", "hsa-miR-1266-5p",
  "hsa-miR-1267", "hsa-miR-1268a", "hsa-miR-1268b",
  "hsa-miR-1269a", "hsa-miR-1269b", "hsa-miR-1270",
  "hsa-miR-1271-5p", "hsa-miR-1272", "hsa-miR-1273a",
  "hsa-miR-1273c", "hsa-miR-1273e", "hsa-miR-1273f",
  "hsa-miR-1273g-3p", "hsa-miR-1275", "hsa-miR-1276",
  "hsa-miR-1277-5p", "hsa-miR-1278", "hsa-miR-1279",
  "hsa-miR-1280", "hsa-miR-1281", "hsa-miR-1282",
  "hsa-miR-1283", "hsa-miR-1284", "hsa-miR-1285-5p",
  "hsa-miR-1286", "hsa-miR-1287-5p", "hsa-miR-1288-3p",
  "hsa-miR-1289", "hsa-miR-1290", "hsa-miR-1291",
  "hsa-miR-1292-3p", "hsa-miR-1293", "hsa-miR-1294",
  "hsa-miR-1295a", "hsa-miR-1295b-5p", "hsa-miR-1296-5p",
  "hsa-miR-1297", "hsa-miR-1298-5p", "hsa-miR-1299",
  "hsa-miR-1300", "hsa-miR-1301-3p", "hsa-miR-1302",
  "hsa-miR-1303", "hsa-miR-1304-5p", "hsa-miR-1305",
  "hsa-miR-1306-5p", "hsa-miR-1307-3p", "hsa-miR-1308",
  "hsa-miR-1468-5p", "hsa-miR-1469", "hsa-miR-1470",
  "hsa-miR-1471", "hsa-miR-1537-3p", "hsa-miR-1538",
  "hsa-miR-1539", "hsa-miR-1825", "hsa-miR-1826",
  "hsa-miR-1827", "hsa-miR-1908-5p", "hsa-miR-1909-3p",
  "hsa-miR-1910-5p", "hsa-miR-1911-5p", "hsa-miR-1912",
  "hsa-miR-1913", "hsa-miR-1914-5p", "hsa-miR-1915-3p",
  "hsa-miR-1972", "hsa-miR-1973", "hsa-miR-1976",
  "hsa-miR-2110", "hsa-miR-2113", "hsa-miR-2114-5p",
  "hsa-miR-2276", "hsa-miR-2277-3p", "hsa-miR-2278",
  "hsa-miR-2355-5p", "hsa-miR-2467-3p", "hsa-miR-2681",
  "hsa-miR-2682-5p", "hsa-miR-2861", "hsa-miR-2909",
  "hsa-miR-3065-5p", "hsa-miR-3074-5p", "hsa-miR-3115",
  "hsa-miR-3120-3p", "hsa-miR-3125", "hsa-miR-3126-5p",
  "hsa-miR-3127-5p", "hsa-miR-3128", "hsa-miR-3129-5p",
  "hsa-miR-3130-5p", "hsa-miR-3131", "hsa-miR-3132",
  "hsa-miR-3133", "hsa-miR-3134", "hsa-miR-3135a",
  "hsa-miR-3135b", "hsa-miR-3136-5p", "hsa-miR-3137",
  "hsa-miR-3138", "hsa-miR-3139", "hsa-miR-3140-3p",
  "hsa-miR-3141", "hsa-miR-3142", "hsa-miR-3143",
  "hsa-miR-3144-3p", "hsa-miR-3145-3p", "hsa-miR-3146",
  "hsa-miR-3147", "hsa-miR-3148", "hsa-miR-3149",
  "hsa-miR-3150a-3p", "hsa-miR-3150b-3p", "hsa-miR-3151-5p",
  "hsa-miR-3152-5p", "hsa-miR-3153", "hsa-miR-3154",
  "hsa-miR-3155a", "hsa-miR-3155b", "hsa-miR-3156-5p",
  "hsa-miR-3157-5p", "hsa-miR-3158-5p", "hsa-miR-3159",
  "hsa-miR-3160-5p", "hsa-miR-3161", "hsa-miR-3162-5p",
  "hsa-miR-3163", "hsa-miR-3164", "hsa-miR-3165",
  "hsa-miR-3166", "hsa-miR-3167", "hsa-miR-3168",
  "hsa-miR-3169", "hsa-miR-3170", "hsa-miR-3171",
  "hsa-miR-3173-5p", "hsa-miR-3174", "hsa-miR-3175",
  "hsa-miR-3176", "hsa-miR-3177-3p", "hsa-miR-3178",
  "hsa-miR-3179", "hsa-miR-3180-3p", "hsa-miR-3181",
  "hsa-miR-3182", "hsa-miR-3183", "hsa-miR-3184-5p",
  "hsa-miR-3185", "hsa-miR-3186-3p", "hsa-miR-3187-5p",
  "hsa-miR-3188", "hsa-miR-3189-3p", "hsa-miR-3190-5p",
  "hsa-miR-3191-3p", "hsa-miR-3192", "hsa-miR-3193",
  "hsa-miR-3194-5p", "hsa-miR-3195", "hsa-miR-3196",
  "hsa-miR-3197", "hsa-miR-3198", "hsa-miR-3199",
  "hsa-miR-3200-5p", "hsa-miR-3201", "hsa-miR-3202",
  "hsa-miR-3529-3p", "hsa-miR-3613-5p", "hsa-miR-3614-5p",
  "hsa-miR-3615", "hsa-miR-3616-5p", "hsa-miR-3617",
  "hsa-miR-3618", "hsa-miR-3619-5p", "hsa-miR-3620-5p",
  "hsa-miR-3621", "hsa-miR-3622a-3p", "hsa-miR-3622b-3p",
  "hsa-miR-3646", "hsa-miR-3648", "hsa-miR-3649",
  "hsa-miR-3650", "hsa-miR-3651", "hsa-miR-3652",
  "hsa-miR-3653", "hsa-miR-3654", "hsa-miR-3655",
  "hsa-miR-3656", "hsa-miR-3660", "hsa-miR-3661",
  "hsa-miR-3662", "hsa-miR-3663-3p", "hsa-miR-3664-3p",
  "hsa-miR-3665", "hsa-miR-3666", "hsa-miR-3667-5p",
  "hsa-miR-3668", "hsa-miR-3669", "hsa-miR-3670",
  "hsa-miR-3671", "hsa-miR-3673", "hsa-miR-3674",
  "hsa-miR-3675", "hsa-miR-3676", "hsa-miR-3677-3p",
  "hsa-miR-3678-3p", "hsa-miR-3679-5p", "hsa-miR-3680-3p",
  "hsa-miR-3681-5p", "hsa-miR-3682-3p", "hsa-miR-3683",
  "hsa-miR-3684", "hsa-miR-3685", "hsa-miR-3686",
  "hsa-miR-3687", "hsa-miR-3688-3p", "hsa-miR-3689a-5p",
  "hsa-miR-3689b-5p", "hsa-miR-3689c", "hsa-miR-3689d-5p",
  "hsa-miR-3689e", "hsa-miR-3690-3p", "hsa-miR-3691-5p",
  "hsa-miR-3692-5p", "hsa-miR-3907", "hsa-miR-3908",
  "hsa-miR-3909", "hsa-miR-3910", "hsa-miR-3911",
  "hsa-miR-3912-3p", "hsa-miR-3913-5p", "hsa-miR-3914",
  "hsa-miR-3915", "hsa-miR-3916", "hsa-miR-3917",
  "hsa-miR-3918", "hsa-miR-3919", "hsa-miR-3920",
  "hsa-miR-3921", "hsa-miR-3922-3p", "hsa-miR-3923",
  "hsa-miR-3924", "hsa-miR-3925-5p", "hsa-miR-3926-3p",
  "hsa-miR-3927-3p", "hsa-miR-3928-3p", "hsa-miR-3929",
  "hsa-miR-3931", "hsa-miR-3933", "hsa-miR-3934-5p",
  "hsa-miR-3935", "hsa-miR-3936", "hsa-miR-3937",
  "hsa-miR-3938", "hsa-miR-3939", "hsa-miR-3940-5p",
  "hsa-miR-3941", "hsa-miR-3942-5p", "hsa-miR-3943",
  "hsa-miR-3944-5p", "hsa-miR-3945"
)

cat("主要人类miRNA总数:", length(major_mirnas), "\n")

# 分批次查询 (multiMiR 限制每批查询量)
batch_size <- 50
all_results <- list()
n_batches <- ceiling(length(major_mirnas) / batch_size)
success_count <- 0

cat("开始分", n_batches, "批查询 multiMiR...\n")

for (i in 1:n_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(major_mirnas))
  batch <- major_mirnas[start_idx:end_idx]
  
  cat(sprintf("  批次%3d/%d: miRNA %4d-%4d ...", i, n_batches, start_idx, end_idx))
  flush.console()
  
  tryCatch({
    r <- get_multimir(
      mirna = batch,
      org = "hsa",
      table = "validated",
      predicted.cutoff.type = "n",
      use.tibble = FALSE,
      legacy.out = FALSE
    )
    if (!is.null(r@data) && nrow(r@data) > 0) {
      all_results[[i]] <- r@data
      success_count <- success_count + 1
      cat(" OK (", nrow(r@data), "条)\n", sep = "")
    } else {
      cat(" 空结果\n")
    }
  }, error = function(e) {
    cat(" ERROR:", conditionMessage(e), "\n")
  })
  
  Sys.sleep(2)
}

cat("\n成功:", success_count, "/", n_batches, " 批次\n")

if (length(all_results) == 0) {
  cat("ERROR: 未获取到任何数据\n")
  cat("尝试使用备选策略：通过靶基因反向查询...\n")
  quit(status = 1)
}

# 合并
combined <- do.call(rbind, all_results)
cat("合并后总行数:", nrow(combined), "\n")

# 保存所有结果
output_file <- file.path(output_dir, "gene_mirna_edges.txt")

# 检测列名
cols <- colnames(combined)
cat("\n可用列:", paste(cols, collapse = ", "), "\n")

target_col <- grep("target_symbol|target_gene|Target_Gene", cols, value = TRUE)
mirna_col <- grep("mirna|mature_mirna_id", cols, value = TRUE)
experiment_col <- grep("experiment|Support_Type", cols, value = TRUE)

target_col <- target_col[1] %||% "target_symbol"
mirna_col <- mirna_col[1] %||% "mature_mirna_id"
experiment_col <- experiment_col[1] %||% NULL

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a

# 检测实际列名
target_col <- if ("target_symbol" %in% cols) "target_symbol" else
              if ("Target.Gene" %in% cols) "Target.Gene" else
              if ("target_gene" %in% cols) "target_gene" else cols[1]

mirna_col <- if ("mature_mirna_id" %in% cols) "mature_mirna_id" else
             if ("miRNA" %in% cols) "miRNA" else
             if ("mirna_id" %in% cols) "mirna_id" else cols[2]

exp_col <- if ("experiment" %in% cols) "experiment" else
           if ("Support_Type" %in% cols) "Support_Type" else NULL

cat("写入: gene_col='", target_col, "', mirna_col='", mirna_col, "'\n", sep = "")
if (!is.null(exp_col)) cat("  实验证据列: '", exp_col, "'\n", sep = "")

write_cols <- c(target_col, mirna_col)
if (!is.null(exp_col) && exp_col %in% cols) write_cols <- c(write_cols, exp_col)

write.table(
  combined[, write_cols],
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("OK! 保存到:", output_file, "\n")
cat("总边数:", nrow(combined), "\n")
cat("唯一miRNA:", length(unique(combined[[mirna_col]])), "\n")
cat("唯一靶基因:", length(unique(combined[[target_col]])), "\n")

# 保存统计
stats_file <- file.path(output_dir, "gene_mirna_stats.json")
cat("统计文件:", stats_file, "\n")

library(jsonlite)
stats <- list(
  total_edges = nrow(combined),
  unique_mirnas = length(unique(combined[[mirna_col]])),
  unique_genes = length(unique(combined[[target_col]])),
  source = "multiMiR v2.4.0 (miRTarBase validated data)",
  url = "http://multimir.org/"
)
writeLines(toJSON(stats, pretty = TRUE, auto_unbox = TRUE), stats_file)

cat("\n===== 完成! =====\n")