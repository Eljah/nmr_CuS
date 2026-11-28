package org.eljah.nmr.relaxation;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

/** Small dependency-free UTF-8 CSV/JSON utilities; no computation is delegated to scripts. */
final class FilesIO {
    private FilesIO() {}
    static List<Map<String,String>> csv(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.isEmpty()) return List.of();
        List<String> header = fields(lines.get(0).replace("\uFEFF", ""));
        List<Map<String,String>> rows = new ArrayList<>();
        for (int k=1; k<lines.size(); k++) {
            if (lines.get(k).isBlank()) continue;
            List<String> values = fields(lines.get(k));
            Map<String,String> row = new LinkedHashMap<>();
            for (int j=0; j<header.size(); j++)
                row.put(header.get(j), j<values.size()?values.get(j):"");
            rows.add(row);
        }
        return rows;
    }
    static List<String> fields(String line) {
        List<String> result = new ArrayList<>();
        StringBuilder field = new StringBuilder();
        boolean quoted = false;
        for (int j=0; j<line.length(); j++) {
            char c = line.charAt(j);
            if (c=='"') {
                if (quoted && j+1<line.length() && line.charAt(j+1)=='"') { field.append('"'); j++; }
                else quoted=!quoted;
            } else if (c==',' && !quoted) { result.add(field.toString()); field.setLength(0); }
            else field.append(c);
        }
        if (quoted) throw new IllegalArgumentException("Unclosed CSV quote: "+line);
        result.add(field.toString());
        return result;
    }
    static double number(Map<String,String> row, String key) {
        String text=row.getOrDefault(key,"").trim();
        if (text.isEmpty()) return Double.NaN;
        try { return Double.parseDouble(text); }
        catch (NumberFormatException ex) { return Double.NaN; }
    }
    static void writeCsv(Path path, List<String> header, List<List<Object>> rows) throws IOException {
        Files.createDirectories(path.getParent());
        StringBuilder out=new StringBuilder();
        appendRow(out,new ArrayList<>(header));
        for (List<Object> row: rows) appendRow(out,row);
        Files.writeString(path,out,StandardCharsets.UTF_8);
    }
    private static void appendRow(StringBuilder out,List<?> values) {
        for (int j=0;j<values.size();j++) {
            if(j>0)out.append(',');
            Object value=values.get(j);
            String s=value==null?"":value.toString();
            if(value instanceof Double d && !Double.isFinite(d))s="";
            if(s.indexOf(',')>=0 || s.indexOf('"')>=0 || s.indexOf('\n')>=0)
                out.append('"').append(s.replace("\"","\"\"")).append('"');
            else out.append(s);
        }
        out.append('\n');
    }
    static List<Object> row(Object...values){return new ArrayList<>(Arrays.asList(values));}
    static Map<String,Object> map(Object...pairs) {
        Map<String,Object> result=new LinkedHashMap<>();
        for(int j=0;j<pairs.length;j+=2)result.put(pairs[j].toString(),pairs[j+1]);
        return result;
    }
    static void writeJson(Path path,Object object)throws IOException {
        Files.createDirectories(path.getParent());
        Files.writeString(path,json(object,0)+"\n",StandardCharsets.UTF_8);
    }
    static String json(Object o,int depth) {
        if(o==null)return "null";
        if(o instanceof Number n) {
            if(n instanceof Double && !Double.isFinite(n.doubleValue()))return "null";
            return n.toString();
        }
        if(o instanceof Boolean)return o.toString();
        if(o instanceof Map<?,?> m) {
            List<String> parts=new ArrayList<>();
            for(var e:m.entrySet())parts.add("  ".repeat(depth+1)+quote(e.getKey().toString())+": "+json(e.getValue(),depth+1));
            return "{\n"+String.join(",\n",parts)+"\n"+"  ".repeat(depth)+"}";
        }
        if(o instanceof Collection<?> c) {
            List<String> parts=new ArrayList<>();
            for(Object v:c)parts.add(json(v,depth+1));
            return "["+String.join(", ",parts)+"]";
        }
        if(o instanceof double[] a) {
            List<Double> values=new ArrayList<>();for(double v:a)values.add(v);return json(values,depth);
        }
        return quote(o.toString());
    }
    private static String quote(String s) {
        return "\""+s.replace("\\","\\\\").replace("\"","\\\"").replace("\n","\\n")
                .replace("\r","\\r").replace("\t","\\t")+"\"";
    }
    static String sha256(byte[] bytes) {
        try{return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(bytes));}
        catch(NoSuchAlgorithmException e){throw new IllegalStateException(e);}
    }
}
