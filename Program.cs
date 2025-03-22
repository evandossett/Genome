using System.Collections.Concurrent;
using System.Configuration;
using System.Text.RegularExpressions;

internal class Program {
	private static readonly string[] chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "x", "y", "MT"];
	private static readonly List<SNP>[] chromosomalSNPS = new List<SNP>[chromosomes.Length];
	private static readonly string[][] chromosomalRaw = new string[chromosomes.Length][];
	private static readonly db_Chromosome[] db_chromosomes = new db_Chromosome[chromosomes.Length];
	private static string filePath = Environment.CurrentDirectory + "\\";
	private static string genomeFileName = ConfigurationManager.AppSettings["genomeFileName"] ?? "genome.txt";
	private static readonly int maxRetries = 10;
	private static int curRetries = 0;
	private static readonly bool useSigFilter = bool.Parse(ConfigurationManager.AppSettings["useSigFilter"] ?? "true");
	private static readonly string sigFilterExtra = ConfigurationManager.AppSettings["sigFilterExtra"] ?? "benign,not-provided";
	private static readonly int requestsPerSecond = 3;
	private static int cur_requestsPerSecond = 0;
	private static readonly bool reuseDbFiles = bool.Parse(ConfigurationManager.AppSettings["reuseDbFiles"] ?? "true");
	private static readonly bool checkForUpdates = bool.Parse(ConfigurationManager.AppSettings["checkForUpdates"] ?? "true");
	private static readonly int returnMax = int.Parse(ConfigurationManager.AppSettings["checkForUpdates"] ?? "1000");
	private static int done = 0;

	private static void Main(string[] args) {
		string[] raw_genome = [];
		try {
			raw_genome = File.ReadAllLinesAsync(filePath + genomeFileName).Result;

			if (raw_genome == Array.Empty<string>())
				throw new Exception("Genome file was empty");
		}
		catch { Console.WriteLine("Failed reading genome file, ensure that it is in the application folder, matches the file name entered into the App.config, and is in tvs {snpID  chromosome  allele}"); throw; }

		List<SNP> SNPs = RawGenomeToSNP(raw_genome);
		MatchSNPstoChromosomes(SNPs);
		GetDataBase();
		PrepareData();
		SaveData();

		Console.WriteLine("Press enter to exit...");
		Console.ReadLine();
		Environment.Exit(0);
	}

	#region Main Loop
	static List<SNP> RawGenomeToSNP(string[] raw_genome) {
		Console.Write("Converting Raw Data to SNP List: Starting                    ");
		Console.Write("\rConverting Raw Data to SNP List: In Progress                    ");

		ConcurrentBag<SNP> _SNPs = new ConcurrentBag<SNP>();
		Parallel.ForEach(raw_genome, s => {
			if (s.StartsWith('#') || s.StartsWith('i'))
				return;

			string[] split;
			split = s.Split('\t');
			if (split.Length != 4) {
				Console.Write("\rConverting Raw Data to SNP List: Failed                    \n");
				throw new Exception("Failed Parsing on: {" + s + "}");
			}

			Tuple<string, string> id_split = seperateRSIDandPrefix(split[0]);
			SNP snp = new() {
				id = id_split.Item1,
				id_pref = id_split.Item2,
				chromosome = split[1],
				position = split[2],
				genotype = GenotypeToAlleleCharArray(split[3])
			};

			_SNPs.Add(snp);
		});

		Console.Write("\rConverting Raw Data to SNP List: Done                    \n");
		return _SNPs.ToList();
	}

	static void MatchSNPstoChromosomes(List<SNP> SNPs) {
		Console.Write("Sorting SNP List into Respective Chromosone List: Starting                    ");
		Console.Write("\rSorting SNP List into Respective Chromosone List: In Progress                    ");

		for (int i = 0 ; i < chromosomes.Length ; i++) {
			ConcurrentBag<SNP> snp = new ConcurrentBag<SNP>();
			Parallel.ForEach(SNPs, s => {
				if (chromosomes[i].Equals(s.chromosome)) {
					snp.Add(s);
					return;
				}
			});
			chromosomalSNPS[i] = snp.ToList();
		}

		Console.Write("\rSorting SNP List into Respective Chromosone List: Done                    \n");
	}

	private static void GetDataBase() {

		Console.Write("Getting Data from e-utils: Starting                                                                         ");

		List<int> unfinishedChromosomeIndexes = new List<int>();

		if (reuseDbFiles) {
			ConcurrentBag<int> unfinishedChromosomes = new ConcurrentBag<int>();

			Console.Write("\rGetting Data from e-utils: Finding existing data                                                                           ");

			Parallel.For(0, chromosomes.Length, chromIndex => {
				List<string> finishedSplits = new List<string>();
				if (File.Exists(filePath + $"results\\{chromosomes[chromIndex]}_raw_search_results.xml")) {
					int count = int.Parse(Regex.Match(File.ReadAllText(filePath + $"results\\{chromosomes[chromIndex]}_raw_search_results.xml"), "<Count>(\\d+)</Count>").Groups[1].Value);
					for (int splitIndex = 0 ; splitIndex < count ; splitIndex += returnMax) {
						if (!File.Exists(filePath + $"results\\{chromosomes[chromIndex]}_{splitIndex / returnMax}_raw_summarize_results.xml")) {
							unfinishedChromosomes.Add(chromIndex);
							return;
						}
						else {
							finishedSplits.Add(File.ReadAllText(filePath + $"results\\{chromosomes[chromIndex]}_{splitIndex / returnMax}_raw_summarize_results.xml"));
						}

					}
				}
				else {
					unfinishedChromosomes.Add(chromIndex);
					return;
				}

				chromosomalRaw[chromIndex] = finishedSplits.ToArray();
			});

			unfinishedChromosomeIndexes.AddRange(unfinishedChromosomes);
			unfinishedChromosomeIndexes = unfinishedChromosomeIndexes.OrderBy(i => i).ToList();

			string allUnfinished = "";
			foreach (int i in unfinishedChromosomeIndexes) {
				allUnfinished += $"C{i}, ";
			}
			if (unfinishedChromosomeIndexes.Count > 0) {
				allUnfinished = "\b\b";
			}

			Console.Write($"\rGetting Data from e-utils: Did not find files for {allUnfinished}                                ");
		}
		else {
			for (int i = 0 ; i < chromosomes.Length ; i++)
				unfinishedChromosomeIndexes.Add(i);
		}

		Console.Write("\rGetting Data from e-utils: In Progress                                                              ");

		for (int chromosomeIndex = 0 ; chromosomeIndex < chromosomes.Length ; chromosomeIndex++) {
			if (unfinishedChromosomeIndexes.Contains(chromosomeIndex)) {
				Console.Write("\rGetting Data from e-utils: Downloading (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
				chromosomalRaw[chromosomeIndex] = GetChromosomeData(chromosomeIndex);
				Console.Write("\rGetting Data from e-utils: Downloaded (Chromosome: " + chromosomes[chromosomeIndex] + ")                         ");
			}
			else if (checkForUpdates) {
				Console.Write("\rGetting Data from e-utils: Checking for update (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
				if (NeedsUpdate(chromosomeIndex)) {
					Console.Write("\rGetting Data from e-utils: Update Needed (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
					Console.Write("\rGetting Data from e-utils: Downloading Update (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
					chromosomalRaw[chromosomeIndex] = GetChromosomeData(chromosomeIndex);
					Console.Write("\rGetting Data from e-utils: Downloaded Update (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
				}
				else {
					Console.Write("\rGetting Data from e-utils: No Update Needed (Chromosome: " + chromosomes[chromosomeIndex] + ")                        ");
				}
			}
		}

		Console.Write("\rGetting Data from e-utils: Done                            \n");
	}

	private static void PrepareData() {
		Console.Write("Preparing Data: Starting                                ");
		Console.Write("\rPreparing Data: In Progress                           ");
		Console.Write("\rPreparing Data: Deserializing XML Data                ");
		done = 0;
		Parallel.For(0, chromosomes.Length, (i) => { DeserializeDatabaseXML(i); });
		Console.Write("\rPreparing Data: Finished Deserializing                ");
		Console.Write("\rPreparing Data: Filtering Data                        ");
		done = 0;
		Parallel.For(0, chromosomes.Length, (i) => { FilterData(i); });
		Console.Write("\rPreparing Data: Done                                  \n");
	}

	private static void SaveData() {
		Console.Write("Saving Data to Drive: Starting                                           ");
		Console.Write("\rSaving Data to Drive: In Progress                                          ");

		for (int i = 0 ; i < chromosomes.Length ; i++) {
			Console.Write($"\rSaving Data to Drive: Cleaning Previous Data for Chromosome {db_chromosomes[i].id}             ");
			File.Delete(filePath + $"results\\{db_chromosomes[i].id}_filtered_results.txt");
			Console.Write($"\rSaving Data to Drive: Finished Cleaning Previous Data for Chromosome {db_chromosomes[i].id}             ");
		}

		for (int i = 0 ; i < db_chromosomes.Length ; i++) {
			Console.Write($"\rSaving Data to Drive: Begining Save for Chromosome {db_chromosomes[i].id}                     ");
			File.WriteAllText(filePath + $"results\\{db_chromosomes[i].id}_filtered_results.txt", SerializeChromosome(db_chromosomes[i]));
			Console.Write($"\rSaving Data to Drive: Finished Save for Chromosome {db_chromosomes[i].id}                   ");
			db_chromosomes[i] = new db_Chromosome();
		}

		Console.Write("\rSaving Data to Drive: Done                                     \n");
	}
	#endregion

	#region Raw Genome to SNP
	static char[] GenotypeToAlleleCharArray(string Genotype) {
		if (Genotype.IndexOf('-') == -1)
			return Genotype.ToCharArray();
		return Genotype.Remove(0, Genotype.IndexOf('-')).ToCharArray();
	}

	static Tuple<string, string> seperateRSIDandPrefix(string RSID) {
		if (Char.IsNumber(RSID[0]))
			return Tuple.Create(RSID, "");

		for (int i = 1 ; i < RSID.Length ; i++) {
			if (Char.IsNumber(RSID[i])) {
				return Tuple.Create(RSID.Remove(0, i), RSID.Remove(i));
			}
		}

		throw new Exception("RSID Parsing Fail");
	}
	#endregion

	#region Get Database
	private static string[] GetChromosomeData(int ChromID) {
		Directory.CreateDirectory(filePath + $"\\{chromosomes[ChromID]}");
		Console.Write($"\rGetting Data from e-utils: Assembling URL (Chromosome: {chromosomes[ChromID]})                     ");
		string baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
		string searchUrl = $"{baseUrl}esearch.fcgi?";
		string esummaryUrl = $"{baseUrl}efetch.fcgi?";
		string result = "";
		List<string> results = new List<string>();
		using (HttpClient client = new HttpClient()) {
			Console.Write($"\rGetting Data from e-utils: Splitting Filters (Chromosome: {chromosomes[ChromID]})                      ");
			Console.Write($"\rGetting Data from e-utils: Searching (Chromosome: {chromosomes[ChromID]})                      ");
			result = AttemptSearchTask(searchUrl + $"db=snp&term=\"chr%20{chromosomes[ChromID]}%20snp\"%5BFilter%5D+AND+\"snp%20clinvar\"%5BFilter%5D&usehistory=y", client);
			File.Delete(filePath + $"results\\{chromosomes[ChromID]}_raw_seach_results.xml");
			File.WriteAllText(filePath + $"results\\{chromosomes[ChromID]}_raw_search_results.xml", result);
			string webEnv = Regex.Match(result, "<WebEnv>(\\S+)</WebEnv>").Groups[1].Value;
			string queryKey = Regex.Match(result, "<QueryKey>(\\d+)</QueryKey>").Groups[1].Value;
			int count = int.Parse(Regex.Match(result, "<Count>(\\d+)</Count>").Groups[1].Value);

			for (int i = 0 ; i < count ; i += returnMax) {
				if (File.Exists(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml") && reuseDbFiles) {
					results.Add(File.ReadAllText(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml"));
				}
				else {
					esummaryUrl = $"{baseUrl}efetch.fcgi?db=snp&WebEnv={webEnv}&query_key={queryKey}&retstart={i}&retmax={returnMax}&rettype=fasta&retmode=text";
					result = AttemptSummarizeTask(esummaryUrl, client);
					results.Add(result);
					File.Delete(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml");
					File.WriteAllText(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml", result);
				}
			}
		}

		return results.ToArray();
	}

	private static bool NeedsUpdate(int ChromID) {
		Directory.CreateDirectory(filePath + $"\\{chromosomes[ChromID]}");
		Console.Write($"\rGetting Data from e-utils: Assembling URL (Chromosome: {chromosomes[ChromID]})                     ");
		string baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
		string searchUrl = $"{baseUrl}esearch.fcgi?";
		string esummaryUrl = $"{baseUrl}efetch.fcgi?";
		string result = "";
		List<string> results = new List<string>();
		using (HttpClient client = new HttpClient()) {
			Console.Write($"\rGetting Data from e-utils: Splitting Filters (Chromosome: {chromosomes[ChromID]})                      ");
			Console.Write($"\rGetting Data from e-utils: Searching (Chromosome: {chromosomes[ChromID]})                      ");
			result = AttemptSearchTask(searchUrl + $"db=snp&term=\"chr%20{chromosomes[ChromID]}%20snp\"%5BFilter%5D+AND+\"snp%20clinvar\"%5BFilter%5D&usehistory=y", client);
			File.Delete(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_seach_results.xml");
			File.WriteAllText(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_search_results.xml", result);
			string webEnv = Regex.Match(result, "<WebEnv>(\\S+)</WebEnv>").Groups[1].Value;
			string queryKey = Regex.Match(result, "<QueryKey>(\\d+)</QueryKey>").Groups[1].Value;
			int count = int.Parse(Regex.Match(result, "<Count>(\\d+)</Count>").Groups[1].Value);

			if (File.Exists(filePath + $"results\\{chromosomes[ChromID]}_raw_search_results.xml")) {
				int countComparison = int.Parse(Regex.Match(File.ReadAllText(filePath + $"results\\{chromosomes[ChromID]}_raw_search_results.xml"), "<Count>(\\d+)</Count>").Groups[1].Value);
				if (countComparison == count) {
					return false;
				}
				else {
					File.Delete(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_seach_results.xml");
					File.Delete(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_seach_results.xml");
					for (int i = 0 ; File.Exists(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml") ; i++) {
						File.Delete(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml");
					}
					return true;
				}
			}
			else {
				File.Delete(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_seach_results.xml");
				File.Delete(filePath + $"results\\{chromosomes[ChromID]}_comparison_raw_seach_results.xml");
				for (int i = 0 ; File.Exists(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml") ; i++) {
					File.Delete(filePath + $"results\\{chromosomes[ChromID]}_{i / returnMax}_raw_summarize_results.xml");
				}
				return true;
			}
		}
	}

	private static string AttemptSearchTask(string URL, HttpClient client) {
		Task<string> task;
		try {
			if (cur_requestsPerSecond >= requestsPerSecond) {
				Thread.Sleep(1000);
				cur_requestsPerSecond = 0;
			}
			Console.Write($"\rGetting Data from e-utils: Posting Content                    ");
			cur_requestsPerSecond++;
			task = client.GetStringAsync(URL);
			Console.Write($"\rGetting Data from e-utils: Waiting For Search Response                    ");
			task.Wait();
			Console.Write($"\rGetting Data from e-utils: Server Responded                    ");
			curRetries = 0;
			return task.Result;
		}
		catch (Exception ex) {
			if (!ex.Message.Contains("429") && !ex.Message.Contains("408") && !ex.Message.Contains("403"))
				throw;
			else if (curRetries >= maxRetries)
				throw new Exception("Reached Max Retry Limit on GetTask");

			curRetries++;
			return AttemptSearchTask(URL, client);
		}
	}

	private static string AttemptSummarizeTask(string URL, HttpClient client) {
		Task<string> task;
		try {
			if (cur_requestsPerSecond >= requestsPerSecond) {
				Thread.Sleep(1000);
				cur_requestsPerSecond = 0;
			}
			Console.Write($"\rGetting Data from e-utils: Getting Summary String                    ");
			cur_requestsPerSecond++;
			task = client.GetStringAsync(URL);
			Console.Write($"\rGetting Data from e-utils: Waiting For Summary Response                    ");
			task.Wait();
			Console.Write($"\rGetting Data from e-utils: Server Responded                    ");
			curRetries = 0;
			return task.Result;
		}
		catch (Exception ex) {
			if (!ex.Message.Contains("429") && !ex.Message.Contains("408") && !ex.Message.Contains("403"))
				throw;
			else if (curRetries >= maxRetries)
				throw new Exception("Reached Max Retry Limit on GetTask");

			curRetries++;
			return AttemptSummarizeTask(URL, client);
		}
	}
	#endregion

	#region Prepare Data
	private static void DeserializeDatabaseXML(int chromosomeIndex) {
		foreach (string raw_data in chromosomalRaw[chromosomeIndex]) {
			string[] snps_raw = Regex.Matches(raw_data, "<DocumentSummary uid=\"(.*?)\">(.*?)</DocumentSummary>").Select(r => r.Value).ToArray();

			ConcurrentBag<db_SNP> list = new ConcurrentBag<db_SNP>();
			Parallel.ForEach(snps_raw, snp => {
				db_SNP db_SNP = new() {
					ID = Regex.Match(snp, "<SNP_ID>(.*?)</SNP_ID>").Groups[1].Value,
					Chromosome = Regex.Match(snp, "<CHR>(.*?)</CHR>").Groups[1].Value,
					Position = Regex.Match(snp, "<CHRPOS_PREV_ASSM>(.*?)</CHRPOS_PREV_ASSM>").Groups[1].Value,
					SNP_Class = Regex.Match(snp, "<SNP_CLASS>(.*?)</SNP_CLASS>").Groups[1].Value
				};

				Match match_alleles = Regex.Match(snp, "<SPDI>(.*?)</SPDI>");
				string[] split = match_alleles.Groups[1].Value.Split([':', ',']);
				List<Tuple<string, string>> alleles = new List<Tuple<string, string>>();
				for (int i = 3 ; i < split.Length ; i += 4) {
					Tuple<string, string> t_alleles = new Tuple<string, string>(split[i - 1].Trim(), split[i].Trim());
					alleles.Add(t_alleles);
				}

				db_SNP.alleles = alleles;

				Match match_MAFs = Regex.Match(snp, "<GLOBAL_MAFS>(.*?)</GLOBAL_MAFS>");
				if (match_MAFs.Groups[1].Success) {
					db_SNP.MAFs = new List<Tuple<string, string>>();
					MatchCollection matches = Regex.Matches(match_MAFs.Groups[1].Value, "<MAF>(.*?)</MAF>");
					foreach (Match m in matches) {
						if (!m.Groups[1].Success)
							continue;

						string Study = "";
						Match MAF_match = Regex.Match(m.Groups[1].Value, "<STUDY>(.*?)</STUDY>");
						if (MAF_match.Groups[1].Success)
							Study = MAF_match.Groups[1].Value;
						string Freq = "";
						MAF_match = Regex.Match(m.Groups[1].Value, "<FREQ>(.*?)</FREQ>");
						if (MAF_match.Success)
							Freq = MAF_match.Groups[1].Value;

						db_SNP.MAFs.Add(Tuple.Create(Study, Freq));
					}
				}

				Match match_Genes = Regex.Match(snp, "<GENES>(.*?)</GENES>");
				if (match_Genes.Groups[1].Success) {
					db_SNP.Genes = new List<Tuple<string, string>>();
					MatchCollection matches = Regex.Matches(match_Genes.Groups[1].Value, "<GENE_E>(.*?)</GENE_E>");
					foreach (Match m in matches) {
						if (!m.Groups[1].Success)
							continue;

						string Name = "";
						Match Genes_match = Regex.Match(m.Groups[1].Value, "<NAME>(.*?)</NAME>");
						if (Genes_match.Success)
							Name = Genes_match.Groups[1].Value;
						string GeneID = "";
						Genes_match = Regex.Match(m.Groups[1].Value, "<GENE_ID>(.*?)</GENE_ID>");
						if (Genes_match.Groups[1].Success)
							GeneID = Genes_match.Groups[1].Value;

						db_SNP.Genes.Add(Tuple.Create(Name, GeneID));
					}
				}

				string Population = "N/A";
				Match match_Population = Regex.Match(snp, "<GLOBAL_POPULATION>(.*?)</GLOBAL_POPULATION>");
				if (match_Population.Groups[1].Success)
					Population = match_Population.Groups[1].Value;
				db_SNP.Population = Population;

				string Clinical_Significance = "N/A";
				Match match_Clinical_Significance = Regex.Match(snp, "<CLINICAL_SIGNIFICANCE>(.*?)</CLINICAL_SIGNIFICANCE>");
				if (match_Clinical_Significance.Groups[1].Success)
					Clinical_Significance = match_Clinical_Significance.Groups[1].Value;
				db_SNP.Clinical_Significance = Clinical_Significance;

				list.Add(db_SNP);
			});

			if (db_chromosomes[chromosomeIndex] == null)
				db_chromosomes[chromosomeIndex] = new db_Chromosome(list.ToList());
			else
				db_chromosomes[chromosomeIndex].db_SNPs.AddRange(list);
		}
		chromosomalRaw[chromosomeIndex] = new string[1];
		done++;
		Console.Write($"\rPreparing Data: Deserializing XML Data ({done}/{chromosomes.Length})              ");
	}

	private static void FilterData(int ChromID) {
		ConcurrentBag<db_SNP> db_snps_filtered = new ConcurrentBag<db_SNP>();
		Parallel.ForEach(chromosomalSNPS[ChromID], snp => {
			foreach (db_SNP db_SNP in db_chromosomes[ChromID].db_SNPs) {
				if (db_SNP.ID.Equals(snp.id) && (string.IsNullOrEmpty(db_SNP.Clinical_Significance) || db_SNP.Clinical_Significance.Split(',').All(s1 => sigFilterExtra.Split(',').Any(s2 => s1.Contains(s2)))) && useSigFilter)
					break;
				else if (db_SNP.ID.Equals(snp.id) && db_SNP.alleles != null && snp.genotype != null) {
					foreach (Tuple<string, string> allele in db_SNP.alleles) {
						if (allele.Item2.Length == 1) {
							if (allele.Item2.Any(a1 => snp.genotype.Any(a2 => a1.Equals(a2)))) {
								db_snps_filtered.Add(db_SNP);
								break;
							}
							else {
								continue;
							}
						}
						else if (allele.Item2.Length == 2) {
							if (snp.genotype.Length == 2 && allele.Equals($"{snp.genotype[0]}{snp.genotype[1]}")) {
								db_snps_filtered.Add(db_SNP);
							}
							else {
								continue;
							}
						}
						else
							continue;
					}
					break;
				}
			}
		});

		db_chromosomes[ChromID] = new db_Chromosome(db_snps_filtered.ToList()) { id = chromosomes[ChromID] };
		chromosomalSNPS[ChromID] = new List<SNP>();
		done++;
		Console.Write($"\rPreparing Data: Filtering Data ({done}/{chromosomes.Length})              ");
	}
	#endregion

	#region Save Data
	private static string SerializeChromosome(db_Chromosome db_chrom) {
		string serialized = "";

		serialized += $"Chromosome: {db_chrom.id}\n";

		foreach (db_SNP snp in db_chrom.db_SNPs) {
			serialized += $"\n\n\n\tID: {snp.ID}";
			serialized += $"\n\t\tPosition: {snp.Position}";
			serialized += $"\n\t\tSNP_Class: {snp.SNP_Class}";
			serialized += $"\n\t\tChromosome: {snp.Chromosome}";
			serialized += $"\n\t\tClinical_Signficance: {snp.Clinical_Significance}";
			serialized += $"\n\t\tAlleles: ";
			if (snp.alleles != null) {
				foreach (Tuple<string, string> a in snp.alleles) {
					serialized += $"\n\t\t\t {a.Item1}>{a.Item2}";
				}
			}
			if (snp.Genes != null) {
				serialized += $"\n\t\tGenes: ";
				foreach (Tuple<string, string> tuple in snp.Genes) {
					serialized += $"\n\t\t\t{tuple.Item1}: {tuple.Item2}";
				}
			}

			if (snp.MAFs != null) {
				serialized += $"\n\t\tHAFs: ";
				foreach (Tuple<string, string> tuple in snp.MAFs) {
					serialized += $"\n\t\t\t{tuple.Item1}: {tuple.Item2}";
				}
			}
			serialized += $"\n\t\tPopulation: {snp.Population}\n";
		}

		return serialized;
	}
	#endregion
}