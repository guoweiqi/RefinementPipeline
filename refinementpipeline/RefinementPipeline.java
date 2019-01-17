package refinementpipeline;

import java.net.*;
import java.util.*;
import java.io.*;

class PipelineFunction {

    /**
     * List of structures that have no Homology Models, but we want to use the
     * experimental structures.
     */
    public String getExperimental = "GJB2";
    /**
     * If true, checks against Github before downloading the DVD structure.
     */
    public boolean github;
    /**
     * Stores gene name for each structure.
     */
    public String gene;
    /**
     * Stores percent disorder calculated by PONDR.
     */
    public String percentDisorder;
    /**
     * UniProt link for each gene.
     */
    public URL uniprotObj;
    /**
     * Indicates whether a Protein Model Portal link exists. Important for
     * downloading structures when there is no ProteinModelPortal link, but
     * there is a SwissModel link.
     */
    public String noProteinModelPortal;
    /**
     * Provides Protein Model Portal link for each gene.
     */
    public String proteinModelPortal;
    /**
     * Stores initial repository links from Protein Model Portal page.
     */
    public List<String> proteinModelLinks = new ArrayList<>();
    /**
     * Stores residue ranges for each structure.
     */
    public List<String> residueList = new ArrayList<>();
    /**
     * Stores sequence identities for each structure.
     */
    public List<Double> seqList = new ArrayList<>();
    /**
     * Stores the direct SwissModel/ModBase download link for each structure.
     */
    public List<String> pdbLinks = new ArrayList<>();
    /**
     * Stores total length for each protein.
     */
    public List<String> proteinLengths = new ArrayList<>();
    /**
     * Stores the type of each structure (e.g. monomer, dimer, hexamer).
     */
    public List<String> structureTypes = new ArrayList<>();

    /**
     * Checks against the DVD Github for structures that have already been
     * refined.
     *
     * @param residueList
     * @param structureTypes
     * @param deleteInt
     * @throws MalformedURLException
     * @throws IOException
     */
    public List<Integer> checkGithub(List<Integer> deleteInt) throws MalformedURLException, IOException {
        /**
         * Creates URL object for Github link.
         */
        URL githubLink = new URL("https://github.com/wtollefson/dvd-structures");
        /**
         * Loops through the different structure types for a single gene.
         */
        for (int i = 0; i < structureTypes.size(); i++) {
            BufferedReader br = new BufferedReader(new InputStreamReader(githubLink.openStream()));
            String webpage = null;
            /**
             * Reads the Github URL.
             */
            while ((webpage = br.readLine()) != null) {
                if (webpage.contains("tree/master/" + gene + structureTypes.get(i))) {
                    URL geneLink = new URL("https://github.com/wtollefson/dvd-structures/tree/master/" + gene + structureTypes.get(i));
                    String webpage1 = null;
                    BufferedReader br2 = new BufferedReader(new InputStreamReader(geneLink.openStream()));
                    while ((webpage1 = br2.readLine()) != null) {
                        /**
                         * If the Github already contains the PDB, add the index
                         * of that PDB to deleteInt.
                         */
                        if (webpage1.contains("master/" + gene + structureTypes.get(i) + "/" + residueList.get(i))) {
                            deleteInt.add(i);
                            break;
                        }
                    }
                    break;
                }
            }
        }
        return deleteInt;
    }

    /**
     * Finds gene name from UniProt website.
     *
     * @param entryID Entry ID for a gene found on UniProt.
     * @return
     * @throws MalformedURLException
     * @throws IOException
     */
    public String getGene(String entryID) throws MalformedURLException, IOException {

        String uniprot = "https://www.uniprot.org/uniprot/" + entryID;
        uniprotObj = new URL(uniprot);

        /**
         * Reads source code from UniProt link.
         */
        BufferedReader br5 = new BufferedReader(new InputStreamReader(uniprotObj.openStream()));

        /**
         * Finds gene name from source code of UniProt website.
         */
        String webpage1 = null;
        while ((webpage1 = br5.readLine()) != null) {
            if (webpage1.contains("- Homo sapiens (Human) -")) {
                break;
            }
        }
        /**
         * Split function splits lines from source code into substrings and
         * places them in arrays--use this to find and store specific strings
         * from source code.
         */
        String firstHalf2 = webpage1.split("lang=\"en\"><head><title>")[1];
        gene = firstHalf2.split(" - ")[0];
        gene = gene.toUpperCase();

        System.out.println(gene);

        return gene;
    }

    /**
     * Calculates the percent disorder of a protein based on its FASTA sequence.
     *
     * @param entryID
     * @return
     * @throws MalformedURLException
     * @throws IOException
     */
    public String getIDP(String entryID) throws MalformedURLException, IOException {
        URL fastaURL = new URL("https://www.uniprot.org/uniprot/" + entryID + ".fasta");
        BufferedReader br = new BufferedReader(new InputStreamReader(fastaURL.openStream()));

        String webpage = null;
        String fasta = null;
        while ((webpage = br.readLine()) != null) {
            fasta += webpage + "\n";
        }
        fasta = fasta.split("null")[1];

        try {
            /**
             * Construct data
             */
            String data = URLEncoder.encode("VLXT", "UTF-8") + "=" + URLEncoder.encode("on", "UTF-8");
            data += "&" + URLEncoder.encode("stats", "UTF-8") + "=" + URLEncoder.encode("on", "UTF-8");
            data += "&" + URLEncoder.encode("Sequence", "UTF-8") + "=" + URLEncoder.encode(fasta, "UTF-8");
            data += "&" + URLEncoder.encode("submit_result", "UTF-8") + "=" + URLEncoder.encode("submit", "UTF-8");

            /**
             * Send data
             */
            URL url = new URL("http://www.pondr.com/cgi-bin/pondr.cgi");
            URLConnection conn = url.openConnection();
            conn.setDoOutput(true);
            OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
            wr.write(data);
            wr.flush();

            /**
             * Get the response
             */
            BufferedReader rd = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            String line;
            while ((line = rd.readLine()) != null) {
                if (line.contains("Overall percent disordered: ")) {
                    String linesplit = line.split(": ")[1];
                    percentDisorder = linesplit.split("\t")[0];
                    break;
                }
            }
            wr.close();
            rd.close();
        } catch (IOException e) {
        }
        return percentDisorder;
    }

    /**
     * Identifies whether a Protein Model Portal link exists. If there is no
     * Protein Model Portal but a SwissModel link exists, then the SwissModel
     * link is stored.
     *
     * @throws MalformedURLException
     * @throws IOException
     */
    public void getNoProteinModelPortal(String entryID) throws MalformedURLException, IOException {

        BufferedReader br5 = new BufferedReader(new InputStreamReader(uniprotObj.openStream()));

        String line2;
        /**
         * While loop tells BufferedReader to read in lines from the source
         * until there is nothing left to read.
         */
        while ((line2 = br5.readLine()) != null) {
            if (line2.contains("ProteinModelPortal")) {
                /*
                String[] proteinModelSplit = line2.split("/uniprot/");
                String firstProteinModel = proteinModelSplit[1];
                String[] proteinModelSplit2 = firstProteinModel.split("\" onclick=");
                proteinModelPortal = proteinModelSplit2[0];
                 */
                proteinModelPortal = entryID;
                break;
            } else if (line2.contains("tertiary") && line2.contains("swissmodel") && !line2.contains("ProteinModelPortal")) {
                String[] smrSplit = line2.split("</td><td><a href=\"");
                String firstSMR = smrSplit[1];
                String[] smrSplit2 = firstSMR.split("\">");
                noProteinModelPortal = smrSplit2[0];
                break;
            }
        }
    }

    /**
     * Gets the initial repository links from the "Query Result" page of the
     * Protein Model Portal website.
     *
     * @throws MalformedURLException
     * @throws IOException
     */
    public void getProteinModelLinks() throws MalformedURLException, IOException {

        URL proteinModelObj = new URL("https://www.proteinmodelportal.org/query/uniprot/" + proteinModelPortal);

        BufferedReader br6 = new BufferedReader(new InputStreamReader(proteinModelObj.openStream()));

        /**
         * Finds model links on the ProteinModelPortal website and adds them to
         * an ArrayList.
         */
        String line3;
        String modbase = null;
        String swissmodel = null;
        while ((line3 = br6.readLine()) != null) {
            if (line3.contains("colorElement") && line3.contains("MODBASE")) {
                String[] modbaseSplit = line3.split("provider=");
                String firstModbase = modbaseSplit[1];
                String[] modbaseSplit2 = firstModbase.split("'\\)\"");
                modbase = modbaseSplit2[0];
                modbase = modbase.replaceAll("amp;", "");
                proteinModelLinks.add(modbase);
            } else if (line3.contains("colorElement") && line3.contains("SWISSMODEL")) {
                String[] swissmodelSplit = line3.split("provider=");
                String firstSwissmodel = swissmodelSplit[1];
                String[] swissmodelSplit2 = firstSwissmodel.split("'\\)\"");
                swissmodel = swissmodelSplit2[0];
                swissmodel = swissmodel.replaceAll("amp;", "");
                proteinModelLinks.add(swissmodel);
            }
        }

        /**
         * Takes all SwissModel links from the list and moves them to a new list
         * for deletion.
         */
        List<String> swissModelLinks = new ArrayList<>();
        for (Iterator<String> iter = proteinModelLinks.listIterator(); iter.hasNext();) {
            String a = iter.next();
            if (a.contains("SWISSMODEL")) {
                swissModelLinks.add(a);
                iter.remove();
            }
        }

        /**
         * Deletes all SwissModel links except for one and adds the one link
         * back to the original list. This eliminates redundancy and saves time
         * (since every SwissModel link redirects to the same site).
         */
        if (swissModelLinks.size() > 0) {
            swissModelLinks.subList(1, swissModelLinks.size()).clear();
            proteinModelLinks.add(swissModelLinks.get(0));
        }
    }

    /**
     * Constructs lists that contain the download link, residue range, sequence
     * identity, structure type (monomer, dimer, trimer, etc), and total protein
     * length for each structure available for a gene. Structures with the same
     * residue range and structure type are sorted based on sequence identity.
     *
     * @throws MalformedURLException
     * @throws IOException
     */
    public void buildLists() throws MalformedURLException, IOException {

        if (noProteinModelPortal != null) {
            proteinModelLinks.add("this is a placeholder for when no Protein Model Portal link is provided but there is a Swissmodel link");
        }

        /**
         * Loops through each model link (usually, one ModBase and one
         * SwissModel).
         */
        for (int x = 0; x < proteinModelLinks.size(); x++) {

            String modelDetails = "https://www.proteinmodelportal.org/?pid=modelDetail&provider=" + proteinModelLinks.get(x);
            URL modelObj = new URL(modelDetails);

            BufferedReader br7 = new BufferedReader(new InputStreamReader(modelObj.openStream()));

            /**
             * Locates source code for the "Model Information" page embedded
             * within the "Model Details" page.
             */
            String line4;
            String detail = null;
            while ((line4 = br7.readLine()) != null) {
                if (line4.contains("Loading Model Detail Information")) {
                    String[] detailSplit = line4.split("pmpuid=");
                    String firstDetail = detailSplit[1];
                    String[] detailSplit2 = firstDetail.split("\", \"");
                    detail = detailSplit2[0];
                    break;
                }
            }

            String modelInfo = "https://www.proteinmodelportal.org/?pid=asyncData&vid=structure&pmpuid=" + detail;
            URL infoObj = new URL(modelInfo);

            /**
             * finds link to ModBase/SwissModel page from the model information
             * source code and stores it as the String variable "pageSource".
             */
            BufferedReader br8 = new BufferedReader(new InputStreamReader(infoObj.openStream()));
            String line5;
            String info = null;
            while ((line5 = br8.readLine()) != null) {
                if (line5.contains("img id")) {
                    String[] infoSplit = line5.split("href=\"");
                    String firstInfo = infoSplit[1];
                    String[] infoSplit2 = firstInfo.split("\" class");
                    info = infoSplit2[0];
                    break;
                }
            }

            String pageSource = info;

            if (noProteinModelPortal != null) {
                pageSource = noProteinModelPortal;
            }

            /**
             * This if statement sorts the links between SwissModel and ModBase.
             */
            if (pageSource.contains("swissmodel")) {

                /**
                 * Fixes any issues with redirects.
                 */
                URL fixRedirect = new URL(pageSource);
                String newUrl = null;

                if (!pageSource.contains("https")) {
                    try {
                        HttpURLConnection conn = (HttpURLConnection) fixRedirect.openConnection();
                        conn.setReadTimeout(5000);
                        conn.addRequestProperty("Accept-Language", "en-US,en;q=0.8");
                        conn.addRequestProperty("User-Agent", "Mozilla");
                        conn.addRequestProperty("Referer", "google.com");

                        boolean redirect = false;

                        int status = conn.getResponseCode();
                        if (status != HttpURLConnection.HTTP_OK) {
                            if (status == HttpURLConnection.HTTP_MOVED_TEMP || status == HttpURLConnection.HTTP_MOVED_PERM || status == HttpURLConnection.HTTP_SEE_OTHER) {
                                redirect = true;
                            }
                        }

                        if (redirect) {

                            newUrl = conn.getHeaderField("Location");

                            String cookies = conn.getHeaderField("Set-Cookie");

                            conn = (HttpURLConnection) new URL(newUrl).openConnection();
                            conn.setRequestProperty("Cookie", cookies);
                            conn.addRequestProperty("Accept-Language", "en-US,en;q=0.8");
                            conn.addRequestProperty("User-Agent", "Mozilla");
                            conn.addRequestProperty("Referer", "google.com");
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                    pageSource = newUrl;
                }

                URL pageObj = new URL(pageSource);

                BufferedReader br2 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                /**
                 * Locates every PDB download link under "Homology Models" on
                 * SwissModel and stores them in pdbLinks.
                 */
                String line1;
                while ((line1 = br2.readLine()) != null) {
                    if (line1.contains(">Homology models<")) {
                        line1 = br2.readLine();
                        while (!(line1 = br2.readLine()).contains("tabindex") && !line1.contains("smrSuggest")) {
                            if (line1.contains("repository")) {
                                String[] linkSplit = line1.split("href=\"");
                                String firstLink = linkSplit[2];
                                if (firstLink.contains(".pdb")) {
                                    String[] linkSplit2 = firstLink.split("\" target=");
                                    String link = linkSplit2[0];
                                    pdbLinks.add(link);
                                } else {
                                    linkSplit = line1.split("href=\"");
                                    firstLink = linkSplit[1];
                                    String[] linkSplit2 = firstLink.split("\" target=");
                                    String link = linkSplit2[0];
                                    pdbLinks.add(link);
                                }
                            }
                        }
                    }
                }

                BufferedReader br4 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                /**
                 * Finds the structure type, residue range, sequence identity,
                 * and full protein length from the SwissModel website and
                 * stores them in structureTypes, residueList, seqList, and
                 * proteinLengths. The index for each PDB link should correspond
                 * to the index for each structure type, residue range, etc.
                 */
                String line;
                String residueString = null;
                Double seqId = null;
                String lengthString = null;
                while ((line = br4.readLine()) != null) {
                    if (line.contains(">Homology models<")) {
                        if (line.contains("monomer")) {
                            structureTypes.add("");
                        } else if (line.contains("homo-2-mer")) {
                            structureTypes.add("-Dimer");
                        } else if (line.contains("homo-3-mer")) {
                            structureTypes.add("-Trimer");
                        } else if (line.contains("homo-4-mer")) {
                            structureTypes.add("-Tetramer");
                        } else if (line.contains("homo-5-mer")) {
                            structureTypes.add("-Pentamer");
                        } else if (line.contains("homo-6-mer")) {
                            structureTypes.add("-Hexamer");
                        } else if (line.contains("homo-12-mer")) {
                            structureTypes.add("-12-mer");
                        } else if (line.contains("mer<")) {
                            structureTypes.add("-Heteromer");
                        }
                        /**
                         * The keyword is "tabindex" if the Homology Models on
                         * SwissModel are followed by Homology Models Built on
                         * Isoform Sequence and the keyword is "smrSuggest" if
                         * there are no structures listed after the Homology
                         * Models.
                         */
                        while (!(line = br4.readLine()).contains("tabindex") && !line.contains("smrSuggest")) {
                            if (line.contains("monomer")) {
                                structureTypes.add("");
                            } else if (line.contains("homo-2-mer")) {
                                structureTypes.add("-Dimer");
                            } else if (line.contains("homo-3-mer")) {
                                structureTypes.add("-Trimer");
                            } else if (line.contains("homo-4-mer")) {
                                structureTypes.add("-Tetramer");
                            } else if (line.contains("homo-5-mer")) {
                                structureTypes.add("-Pentamer");
                            } else if (line.contains("homo-6-mer")) {
                                structureTypes.add("-Hexamer");
                            } else if (line.contains("homo-12-mer")) {
                                structureTypes.add("-12-mer");
                            } else if (line.contains("mer<")) {
                                structureTypes.add("-Heteromer");
                            }
                            if (line.contains("seqLength=\"")) {
                                String[] firstSplit = line.split("seqLength=\"");
                                String secondHalf = firstSplit[1];
                                String[] secondSplit = secondHalf.split("\"");
                                lengthString = secondSplit[0];
                                proteinLengths.add(lengthString);
                            }
                            if (line.contains("display:none\">")) {
                                String[] firstSplit = line.split("display:none\">");
                                String secondHalf = firstSplit[1];
                                String[] secondSplit = secondHalf.split("</span>");
                                residueString = secondSplit[0];
                                residueList.add(residueString);
                            }
                            if (line.contains("Range:")) {
                                String[] seqSplit = line.split("</div></td><td>");
                                String firstSeq = seqSplit[1];
                                String[] seqSplit2 = firstSeq.split("</td><td class");
                                String seqString = seqSplit2[0];
                                seqId = Double.parseDouble(seqString);
                                seqList.add(seqId);
                            }
                        }
                    }
                }

                if (gene.contains(getExperimental)) {
                    System.out.println("Grabbing experimental structures");
                    //grab the link for the experimental structure, add to pdbLinks list
                    //include filler data for the sequence identity
                    BufferedReader br3 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                    /**
                     * Locates every PDB download link under "Homology Models"
                     * on SwissModel and stores them in pdbLinks.
                     */
                    String line2;
                    while ((line2 = br3.readLine()) != null) {
                        if (line2.contains(">Experimental structures<")) {
                            line2 = br3.readLine();
                            while (!(line2 = br3.readLine()).contains("tabindex") && !line2.contains("smrSuggest")) {
                                if (line2.contains("repository")) {
                                    String[] linkSplit = line2.split("href=\"");
                                    String firstLink = linkSplit[1];
                                    if (firstLink.contains(".pdb")) {
                                        String[] linkSplit2 = firstLink.split("\" target=");
                                        String link = linkSplit2[0];
                                        pdbLinks.add(link);
                                    } else {
                                        linkSplit = line2.split("href=\"");
                                        firstLink = linkSplit[1];
                                        String[] linkSplit2 = firstLink.split("\" target=");
                                        String link = linkSplit2[0];
                                        pdbLinks.add(link);
                                    }
                                }
                            }
                        }
                    }

                    BufferedReader br5 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                    /**
                     * Finds the structure type, residue range, sequence
                     * identity, and full protein length from the SwissModel
                     * website and stores them in structureTypes, residueList,
                     * seqList, and proteinLengths. The index for each PDB link
                     * should correspond to the index for each structure type,
                     * residue range, etc.
                     */
                    String line3;
                    String residueString1 = null;
                    Double seqId1 = null;
                    String lengthString1 = null;
                    while ((line3 = br5.readLine()) != null) {
                        if (line3.contains(">Experimental structures<")) {
                            if (line3.contains("monomer")) {
                                structureTypes.add("-expt");
                            } else if (line3.contains("homo-2-mer")) {
                                structureTypes.add("-Dimer-expt");
                            } else if (line3.contains("homo-3-mer")) {
                                structureTypes.add("-Trimer-expt");
                            } else if (line3.contains("homo-4-mer")) {
                                structureTypes.add("-Tetramer-expt");
                            } else if (line3.contains("homo-5-mer")) {
                                structureTypes.add("-Pentamer-expt");
                            } else if (line3.contains("homo-6-mer")) {
                                structureTypes.add("-Hexamer-expt");
                            } else if (line3.contains("homo-12-mer")) {
                                structureTypes.add("-12-mer-expt");
                            } else if (line3.contains("mer<")) {
                                structureTypes.add("-Heteromer-expt");
                            }
                            /**
                             * The keyword is "tabindex" if the Homology Models
                             * on SwissModel are followed by Homology Models
                             * Built on Isoform Sequence and the keyword is
                             * "smrSuggest" if there are no structures listed
                             * after the Homology Models.
                             */
                            while (!(line3 = br5.readLine()).contains("tabindex") && !line3.contains("smrSuggest")) {
                                if (line3.contains("monomer")) {
                                    structureTypes.add("-expt");
                                } else if (line3.contains("homo-2-mer")) {
                                    structureTypes.add("-Dimer-expt");
                                } else if (line3.contains("homo-3-mer")) {
                                    structureTypes.add("-Trimer-expt");
                                } else if (line3.contains("homo-4-mer")) {
                                    structureTypes.add("-Tetramer-expt");
                                } else if (line3.contains("homo-5-mer")) {
                                    structureTypes.add("-Pentamer-expt");
                                } else if (line3.contains("homo-6-mer")) {
                                    structureTypes.add("-Hexamer-expt");
                                } else if (line3.contains("homo-12-mer")) {
                                    structureTypes.add("-12-mer-expt");
                                } else if (line3.contains("mer<")) {
                                    structureTypes.add("-Heteromer-expt");
                                }
                                if (line3.contains("seqLength=\"")) {
                                    String[] firstSplit = line3.split("seqLength=\"");
                                    String secondHalf = firstSplit[1];
                                    String[] secondSplit = secondHalf.split("\"");
                                    lengthString1 = secondSplit[0];
                                    proteinLengths.add(lengthString1);
                                }
                                if (line3.contains("display:none\">")) {
                                    String[] firstSplit = line3.split("display:none\">");
                                    String secondHalf = firstSplit[1];
                                    String[] secondSplit = secondHalf.split("</span>");
                                    residueString1 = secondSplit[0];
                                    residueList.add(residueString1);
                                }
                                if (line3.contains("Range:")) {
                                    seqId1 = 0.0; //Experimental structures don't have a sequence identity, 0.0 is used as a placeholder
                                    seqList.add(seqId1);
                                }
                            }
                        }
                    }
                }

            } else if (pageSource.contains("modbase")) {

                URL fixRedirect = new URL(pageSource);
                String newUrl = null;

                try {
                    HttpURLConnection conn = (HttpURLConnection) fixRedirect.openConnection();
                    conn.setReadTimeout(5000);
                    conn.addRequestProperty("Accept-Language", "en-US,en;q=0.8");
                    conn.addRequestProperty("User-Agent", "Mozilla");
                    conn.addRequestProperty("Referer", "google.com");

                    boolean redirect = false;

                    int status = conn.getResponseCode();
                    if (status != HttpURLConnection.HTTP_OK) {
                        if (status == HttpURLConnection.HTTP_MOVED_TEMP || status == HttpURLConnection.HTTP_MOVED_PERM || status == HttpURLConnection.HTTP_SEE_OTHER) {
                            redirect = true;
                        }
                    }

                    if (redirect) {

                        newUrl = conn.getHeaderField("Location");

                        String cookies = conn.getHeaderField("Set-Cookie");

                        conn = (HttpURLConnection) new URL(newUrl).openConnection();
                        conn.setRequestProperty("Cookie", cookies);
                        conn.addRequestProperty("Accept-Language", "en-US,en;q=0.8");
                        conn.addRequestProperty("User-Agent", "Mozilla");
                        conn.addRequestProperty("Referer", "google.com");
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }

                pageSource = newUrl;
                pdbLinks.add(pageSource);

                URL pageObj = new URL(pageSource);

                BufferedReader br4 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                /**
                 * Finds and stores the residue range and sequence identity from
                 * the ModBase website.
                 */
                String line;
                Double seqId = null;
                String residueRange = null;
                while ((line = br4.readLine()) != null) {
                    if (line.contains("Sequence Length</th> <td>")) {
                        String[] firstSplit = line.split("Sequence Length</th> <td>");
                        String secondHalf = firstSplit[1];
                        String[] secondSplit = secondHalf.split("</td></tr></table>");
                        String lengthString = secondSplit[0];
                        proteinLengths.add(lengthString);
                    }
                    if (line.contains("Sequence Identity</th><td><font color=")) {
                        String[] seqSplit = line.split("Sequence Identity</th><td><font color=");
                        String firstSeq = seqSplit[1];
                        String[] seqSplit2 = firstSeq.split("%</font>");
                        String secondSeq = seqSplit2[0];
                        String[] seqSplit3 = secondSeq.split(">");
                        String seqString = seqSplit3[1];
                        seqId = Double.parseDouble(seqString);
                        seqList.add(seqId);
                        structureTypes.add("");
                    }
                    if (line.contains("considered")) {
                        String[] firstSplit = line.split("Target Region</th><td>");
                        String secondHalf = firstSplit[1];
                        String[] secondSplit = secondHalf.split(" </td></tr><tr><th><a");
                        residueRange = secondSplit[0];
                        residueList.add(residueRange);
                    }
                }
            }
        }

        /**
         * If there are two PDB download links that correspond to the same
         * residue range, then one of the PDB links is deleted based on which
         * PDB has the higher sequence identity. The list deleteInt is a list of
         * integers that correspond to the index of the PDB download link that
         * is to be deleted.
         */
        List<Integer> deleteInt = new ArrayList<>();
        /**
         * Nested for loops check for duplicates within residueList.
         */
        for (int j = 0; j < residueList.size(); j++) {
            for (int k = j + 1; k < residueList.size(); k++) {
                if (residueList.get(j).equals(residueList.get(k)) && structureTypes.get(j).equals(structureTypes.get(k))) {
                    if (seqList.get(j) < seqList.get(k)) {
                        deleteInt.add(j);
                    } else if (seqList.get(j) >= seqList.get(k)) {
                        deleteInt.add(k);
                    }
                }
            }
        }

        /**
         * Determines whether to check against Github.
         */
        if (github == true) {
            deleteInt = checkGithub(deleteInt);
        }

        /**
         * Removes duplicates from deleteInt.
         */
        deleteInt = new ArrayList<>(new HashSet<>(deleteInt));

        /**
         * Sorts deleteInt from largest to smallest (this comes in handy when
         * implementing the upcoming for loop).
         */
        Collections.sort(deleteInt, Collections.reverseOrder());

        /**
         * Deletes elements from seqList, residueList, pdbLinks, proteinLengths,
         * and structureTypes based on integers stored in deleteInt (elements in
         * the first five lists should all correspond based on index).
         */
        for (Integer deleteInt1 : deleteInt) {
            int delete = deleteInt1;
            seqList.remove(delete);
            residueList.remove(delete);
            pdbLinks.remove(delete);
            proteinLengths.remove(delete);
            structureTypes.remove(delete);
        }

    }

    /**
     * Downloads PDB files and organizes them into the appropriate directories.
     *
     * @throws IOException
     */
    public void downloadFiles() throws IOException {
        /**
         * Creates a directory titled "pdbFiles".
         */
        File pdbName = new File("pdbFiles");
        String path = pdbName.getCanonicalPath();
        File pdbDir = new File(path);
        pdbDir.mkdir();

        /**
         * Cycles through each link that is now stored in pdbLinks.
         */
        int i;
        for (i = 0; i < pdbLinks.size(); i++) {

            /**
             * Sorts download links between SwissModel and ModBase.
             */
            if (pdbLinks.get(i).contains("repository")) {

                String url = "https://swissmodel.expasy.org" + pdbLinks.get(i);
                /**
                 * The String "file" lists the directory the PDB file will be
                 * downloaded in, followed by the initial name of the file
                 * ("coords.pdb").
                 */
                String file = path + "/coords.pdb";

                BufferedInputStream bufferedIS = null;
                FileOutputStream fileOS = null;

                /**
                 * The try-catch block handles exceptions.
                 */
                try {
                    /**
                     * Downloads the PDB from the download link.
                     */
                    URL urlObj = new URL(url);
                    bufferedIS = new BufferedInputStream(urlObj.openStream());
                    fileOS = new FileOutputStream(file);

                    int data = bufferedIS.read();
                    while (data != -1) {
                        fileOS.write(data);
                        data = bufferedIS.read();
                    }

                    String one = file.replace("coords.pdb", "");

                    String residueRange = residueList.get(i);

                    String underscore = "_";
                    String pdb = ".pdb";
                    String finalgene = null;

                    /**
                     * Creates directories for gene and residue range (following
                     * the conventional notation).
                     */
                    switch (structureTypes.get(i)) {
                        case "-Heteromer":
                            finalgene = gene + "-Heteromer";
                            break;
                        case "-Dimer":
                            finalgene = gene + "-Dimer";
                            break;
                        case "-Trimer":
                            finalgene = gene + "-Trimer";
                            break;
                        case "-Tetramer":
                            finalgene = gene + "-Tetramer";
                            break;
                        case "-Pentamer":
                            finalgene = gene + "-Pentamer";
                            break;
                        case "-Hexamer":
                            finalgene = gene + "-Hexamer";
                            break;
                        case "-12-mer":
                            finalgene = gene + "-12-mer";
                            break;
                        case "-Heteromer-expt":
                            finalgene = gene + "-Heteromer-expt";
                            break;
                        case "-Dimer-expt":
                            finalgene = gene + "-Dimer-expt";
                            break;
                        case "-Trimer-expt":
                            finalgene = gene + "-Trimer-expt";
                            break;
                        case "-Tetramer-expt":
                            finalgene = gene + "-Tetramer-expt";
                            break;
                        case "-Pentamer-expt":
                            finalgene = gene + "-Pentamer-expt";
                            break;
                        case "-Hexamer-expt":
                            finalgene = gene + "-Hexamer-expt";
                            break;
                        case "-12-mer-expt":
                            finalgene = gene + "-12-mer-expt";
                            break;
                        case "-expt":
                            finalgene = gene + "-expt";
                            break;
                        case "":
                            finalgene = gene;
                            break;
                        default:
                            finalgene = gene;
                            break;
                    }

                    File geneDir = new File(one + finalgene);
                    geneDir.mkdir();
                    File residueDir = new File(one + finalgene + "/" + residueRange);
                    residueDir.mkdir();

                    /**
                     * Renames PDB file and copies it to the new directory.
                     */
                    FileInputStream findFile = new FileInputStream(file);
                    String finalName = one + finalgene + "/" + residueRange + "/" + finalgene + underscore + residueRange + pdb;
                    FileOutputStream rename = new FileOutputStream(finalName);

                    int finalData = findFile.read();
                    while (finalData != -1) {
                        rename.write(finalData);
                        finalData = findFile.read();
                    }

                    /**
                     * Deletes the initial PDB file.
                     */
                    File old = new File(file);
                    old.delete();

                    String properties = "SWISSMODEL," + proteinLengths.get(i) + "," + seqList.get(i) + "," + percentDisorder;
                    File propertiesFile = new File(one + finalgene + "/" + residueRange + "/properties.txt");
                    FileWriter out = new FileWriter(propertiesFile, false);
                    out.write(properties);
                    out.flush();
                    out.close();

                } catch (MalformedURLException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                } finally {
                    try {
                        if (fileOS != null) {
                            fileOS.close();
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    try {
                        if (bufferedIS != null) {
                            bufferedIS.close();
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

            } else if (pdbLinks.get(i).contains("modbase")) {

                URL pageObj = new URL(pdbLinks.get(i));

                BufferedReader br2 = new BufferedReader(new InputStreamReader(pageObj.openStream()));

                /**
                 * Finds the different parts of the ModBase download link.
                 */
                String line1;
                String model = null;
                String query = null;
                String human = null;
                while ((line1 = br2.readLine()) != null) {
                    if (line1.contains("email")) {
                        String[] linkSplit = line1.split("\\|");
                        model = linkSplit[2];
                        String[] linkSplit2 = line1.split("\"queryfile\" value=\"");
                        String firstQuery = linkSplit2[1];
                        String[] linkSplit3 = firstQuery.split("\"  /><input");
                        query = linkSplit3[0];
                        String linkSplit4 = linkSplit[3];
                        String[] linkSplit5 = linkSplit4.split("\"  /><p");
                        human = linkSplit5[0];
                        break;
                    }
                }

                /**
                 * Enters the above variables into the template for ModBase
                 * download links.
                 */
                String url = "https://modbase.compbio.ucsf.edu/modbase-cgi/get_modbase_file.cgi/modbase_model_" + model + ".pdb?queryfile=" + query + "=modelid&filetype=model&id=" + model + "|" + human;
                String file = path + "/coords.pdb";

                BufferedInputStream bufferedIS = null;
                FileOutputStream fileOS = null;

                /**
                 * Downloads, renames, and moves PDB file using try-catch block.
                 */
                try {
                    URL urlObj = new URL(url);
                    bufferedIS = new BufferedInputStream(urlObj.openStream());
                    fileOS = new FileOutputStream(file);

                    int data = bufferedIS.read();
                    while (data != -1) {
                        fileOS.write(data);
                        data = bufferedIS.read();
                    }

                    String one = file.replace("coords.pdb", "");

                    String residueRange = residueList.get(i);

                    String underscore = "_";
                    String pdb = ".pdb";

                    File geneDir = new File(one + gene);
                    geneDir.mkdir();
                    File residueDir = new File(one + gene + "/" + residueRange);
                    residueDir.mkdir();

                    FileInputStream findFile = new FileInputStream(file);
                    String finalName = one + gene + "/" + residueRange + "/" + gene + underscore + residueRange + pdb;
                    FileOutputStream rename = new FileOutputStream(finalName);

                    int finalData = findFile.read();
                    while (finalData != -1) {
                        rename.write(finalData);
                        finalData = findFile.read();
                    }
                    File old = new File(file);
                    old.delete();

                    String properties = "MODBASE," + proteinLengths.get(i) + "," + seqList.get(i) + "," + percentDisorder;
                    File propertiesFile = new File(one + gene + "/" + residueRange + "/properties.txt");
                    FileWriter out = new FileWriter(propertiesFile, false);
                    out.write(properties);
                    out.flush();
                    out.close();

                } catch (MalformedURLException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                } finally {
                    try {
                        if (fileOS != null) {
                            fileOS.close();
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    try {
                        if (bufferedIS != null) {
                            bufferedIS.close();
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }
}

public class RefinementPipeline {

    public static void main(String[] args) throws MalformedURLException, IOException {

        /**
         * Loops through different entries entered as arguments on command line.
         */
        final long startTime = System.currentTimeMillis();
        int counter = 0;
        boolean checkGithub = false;

        /**
         * Reads flags.
         */
        for (String arg : args) {
            if (arg.contains("-")) {
                counter++;
                if (arg.equals("-g")) {
                    checkGithub = true;
                }
            }
        }

        for (int i = counter; i < args.length; i++) {
            PipelineFunction pipelineobj = new PipelineFunction();
            pipelineobj.github = checkGithub;
            pipelineobj.getGene(args[i]);
            pipelineobj.getIDP(args[i]);
            pipelineobj.getNoProteinModelPortal(args[i]);
            pipelineobj.getProteinModelLinks();
            pipelineobj.buildLists();
            if (pipelineobj.pdbLinks.isEmpty()) {
                continue;
            }
            pipelineobj.downloadFiles();
        }
        final long endTime = System.currentTimeMillis();

        System.out.println("Finished");
        System.out.println("Total execution time: " + (endTime - startTime));

    }
}
